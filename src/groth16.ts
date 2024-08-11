import { Fr, PointG1, PointG2, validatePairing } from '@kevincharm/blstoise'
import { randomFr } from '@kevincharm/blstoise/dist/src/utils'
import { addPolys, evalPoly } from './poly'
import { QAP, computeSolutionPoly, computeDivisorPoly, tPoly } from './qap'

/**
 * Combined phase 1 & 2 trusted setup outputs for Groth16
 */
export interface TrustedSetup {
    power: number
    tauG1: PointG1[]
    tauG2: PointG2[]
    alphaG1: PointG1
    betaG1: PointG1
    betaG2: PointG2
    gammaG2: PointG2
    deltaG1: PointG1
    deltaG2: PointG2
    cTauPublicG1: PointG1[]
    cTauPrivateG1: PointG1[]
    tTauInvDeltaG1: PointG1[]
}

/**
 * Perform a circuit-specific trusted setup for Groth16 proofs.
 * @param qap Polynomials representing the quadratic arithmetic programme
 */
export function setup({ l, u, v, w }: QAP): TrustedSetup {
    const degree = u[0].length
    const tau = randomFr() // toxic waste
    const alpha = randomFr() // toxic waste
    const beta = randomFr() // toxic waste
    const gamma = randomFr() // toxic waste
    const delta = randomFr() // toxic waste

    const G1 = PointG1.generator()
    const G2 = PointG2.generator()
    const tauG1 = Array.from({ length: degree }, (_, i) => G1.mul(tau.exp(BigInt(i)).value))
    const tauG2 = Array.from({ length: degree }, (_, i) => G2.mul(tau.exp(BigInt(i)).value))
    const alphaG1 = G1.mul(alpha.value)
    const betaG1 = G1.mul(beta.value)
    const betaG2 = G2.mul(beta.value)

    const gammaG2 = G2.mul(gamma.value)
    const deltaG1 = G1.mul(delta.value)
    const deltaG2 = G2.mul(delta.value)

    // foreach row (polynomial) in QAP.u \i -> beta*u_i(x)
    const betaU = u.map((uPoly) => uPoly.map((x) => x.mul(beta)))
    // foreach row (polynomial) in QAP.v \i -> alpha*v_i(x)
    const alphaV = v.map((vPoly) => vPoly.map((x) => x.mul(alpha)))
    // From Groth16, 3.2 NIZK arguments for quadratic arithmetic programs
    // This is the combined sum term from C and AB (for both public and private inputs)
    // evaluated at tau:
    // C = \sum_{i=0}^{m} \beta u_i(\tau) + \alpha v_i(\tau) + w_i(\tau)
    const cG1: PointG1[] = []
    for (let i = 0; i < w.length; i++) {
        const c = G1.mul(evalPoly(addPolys(addPolys(betaU[i], alphaV[i]), w[i]), tau).value)
        cG1.push(c)
    }

    // Public inputs - scale sum term by 1/gamma
    // (1/gamma) * [C]_1
    const invGamma = gamma.inv()
    const cTauPublicG1 = cG1.slice(0, l + 1).map((p) => p.mul(invGamma.value))

    // Private inputs - scale sum term by 1/delta
    // (1/delta) * [C]_1
    const invDelta = delta.inv()
    const cTauPrivateG1 = cG1.slice(l + 1).map((p) => p.mul(invDelta.value))

    // We need this later to compute h(tau)t(tau) in the proving step.
    // See: https://github.com/nikkolasg/playsnark/blob/c406bfdeefbc60345696bfad9105568e91f04c13/groth16.go#L180
    // [tau^i * t(tau) / delta]_1
    const tx = tPoly(degree)
    const tTau = evalPoly(tx, tau)
    const tTauInvDeltaG1 = tauG1.map((ptau) => ptau.mul(tTau.value).mul(invDelta.value))

    return {
        power: degree,
        tauG1,
        tauG2,
        alphaG1,
        betaG1,
        betaG2,
        gammaG2,
        deltaG1,
        deltaG2,
        cTauPublicG1,
        cTauPrivateG1,
        tTauInvDeltaG1,
    }
}

/**
 * Verification key for Groth16 proofs;
 * which is a subset of the circuit-specific trusted setup
 */
export interface VerificationKey {
    alphaG1: PointG1
    betaG2: PointG2
    gammaG2: PointG2
    deltaG2: PointG2
    cTauPublicG1: PointG1[]
}

/** Jens Groth's succinct noninteractive argument of knowledge */
export type Proof = [PointG1, PointG1, PointG2]

/**
 * Compute a succinct Groth16 proof for a given programme and a witness vector.
 * @param setup Circuit-specific trusted setup
 * @param qap Polynomials representing the quadratic arithmetic programme
 * @param witness Witness vector (i.e. a solution to the programme)
 * @param shouldCheckProof Whether to check that the proof verifies before returning
 * @returns A succinct proof that the witness satisfies the QAP
 */
export function prove(
    setup: TrustedSetup,
    qap: QAP,
    witness: Fr[],
    shouldCheckProof?: boolean,
): { proof: Proof; vk: VerificationKey } {
    const r = randomFr() // toxic waste
    const s = randomFr() // toxic waste

    // We fix the first witness to 1 (used for constants)
    const preparedWitness = [Fr.one(), ...witness]
    // Compute polynomials a_i u_i(x)
    const { au, av, ht } = computeSolutionPoly(preparedWitness, qap)

    // A = [\alpha + \sum_{i=0}^{m} a_i u_i(\tau) + r*delta]_1
    const sumA = au.map(({ value }, i) => setup.tauG1[i].mul(value)).reduce((acc, p) => acc.add(p))
    const A = setup.alphaG1.add(sumA).add(setup.deltaG1.mul(r.value))
    // B = [\beta + \sum_{i=0}^{m} a_i v_i(\tau) + s*delta]_1
    const sumBG1 = av
        .map(({ value }, i) => setup.tauG1[i].mul(value))
        .reduce((acc, p) => acc.add(p))
    const BG1 = setup.betaG1.add(sumBG1).add(setup.deltaG1.mul(s.value))
    // B = [\beta + \sum_{i=0}^{m} a_i v_i(\tau) + s*delta]_2
    const sumBG2 = av
        .map(({ value }, i) => setup.tauG2[i].mul(value))
        .reduce((acc, p) => acc.add(p))
    const BG2 = setup.betaG2.add(sumBG2).add(setup.deltaG2.mul(s.value))

    // We need the +h(\tau)*t(\tau) part of C
    const h = computeDivisorPoly(ht, tPoly(qap.u[0].length))
    const htTauG1 = h
        .map((x, i) => setup.tTauInvDeltaG1[i].mul(x.value))
        .reduce((acc, p) => acc.add(p))
    // Private inputs (+ intermediate constraints) part of C
    // \sum_{i=l+1}^{m} a_i(\beta u_i(\tau) + \alpha v_i(\tau) + w_i(\tau)
    const sumC = preparedWitness
        .slice(qap.l + 1)
        .map(({ value }, i) => setup.cTauPrivateG1[i].mul(value))
        .reduce((acc, p) => acc.add(p))
    // -(r * s * \delta)
    const negRsDelta = setup.deltaG1.mul(r.mul(s).value).neg()
    // Full C term evaluated in G1
    // C = \frac{1}{\delta} (\sum_{i=l+1}^{m} a_i(\beta u_i(\tau) + \alpha v_i(\tau) + w_i(\tau)) + h(\tau)t(\tau))
    //  + As + rB - rs\delta
    const C = sumC.add(htTauG1).add(A.mul(s.value)).add(BG1.mul(r.value)).add(negRsDelta)
    // Proof triple, as described by Groth, 3.2 NIZK arguments for quadratic arithmetic programs
    const proof: [PointG1, PointG1, PointG2] = [A, C, BG2]

    // This is the public inputs part
    // \frac{1}{\gamma} (\sum_{i=0}^{l} a_i(\beta u_i(\tau) + \alpha v_i(\tau) + w_i(\tau)))
    const publicInputsSumG1 = setup.cTauPublicG1
        .map((p, i) => p.mul(preparedWitness[i].value))
        .reduce((acc, p) => acc.add(p))

    if (shouldCheckProof) {
        // Check that proof verifies
        const ps = [A.neg(), setup.alphaG1, publicInputsSumG1, C]
        const qs = [BG2, setup.betaG2, setup.gammaG2, setup.deltaG2]
        if (!validatePairing(ps, qs)) {
            throw new Error('Pairing check failed')
        }
    }

    return {
        proof,
        vk: {
            alphaG1: setup.alphaG1,
            betaG2: setup.betaG2,
            gammaG2: setup.gammaG2,
            deltaG2: setup.deltaG2,
            cTauPublicG1: setup.cTauPublicG1,
        },
    }
}

/**
 * Verify a Groth16 proof for a given public input vector.
 *
 * @param proof Groth16 proof
 * @param publicInputs Public input vector
 * @param vk Circuit-specific trusted setup
 * @returns True if the proof is valid
 */
export function verify(proof: Proof, publicInputs: Fr[], vk: VerificationKey): boolean {
    const [A, C, BG2] = proof
    const [alphaG1, betaG2, gammaG2, deltaG2] = [vk.alphaG1, vk.betaG2, vk.gammaG2, vk.deltaG2]

    // We fix the first public input to 1 (used for constants)
    const preparedPublicInputs = [Fr.one(), ...publicInputs]
    if (preparedPublicInputs.length !== vk.cTauPublicG1.length) {
        throw new Error(
            `Invalid number of public inputs: expected ${vk.cTauPublicG1.length - 1}, got ${publicInputs.length}`,
        )
    }
    // \sum_{i=0}^{l} a_i * [\frac{\beta u_i(\tau) + \alpha v_i(\tau) + w_i(\tau)}{\gamma}]_1
    const publicInputsSumG1 = vk.cTauPublicG1
        .map((p, i) => p.mul(preparedPublicInputs[i].value))
        .reduce((acc, p) => acc.add(p))

    // Now from the Groth16 paper, 3.2 NIZK arguments for quadratic arithmetic programs.
    //
    // Parse proof = ([A]_1, [C]_1, [B]_2)
    //
    // Accept the proof iff
    // [A]_1 * [B]_2
    //      = [\alpha]_1 * [\beta]_2
    //      + (\sum_{i=0}^{l} a_i * [\frac{\beta u_i(\tau) + \alpha v_i(\tau) + w_i(\tau)}{\gamma}]_1) * [\gamma]_2
    //      + [C]_1 * [\delta]_2
    // Turn this into a pairing check:
    // Given asymmetric pairing e(P,Q) :: G1 x G2 -> GT
    // e([-A]_1, [B]_2) * e([\alpha]_1, [\beta]_2) * e([publicInputsSum]_1, [\gamma]_2) * e([C]_1, [\delta]_2) = 1
    const ps = [A.neg(), alphaG1, publicInputsSumG1, C]
    const qs = [BG2, betaG2, gammaG2, deltaG2]
    return validatePairing(ps, qs)
}
