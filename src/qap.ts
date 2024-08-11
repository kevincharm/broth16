import { Fr } from '@kevincharm/blstoise'
import { transpose, interpolate, addPolys, mulPolys, subPolys, evalPoly, divPolys } from './poly'

export type FrMatrix = Fr[][]

/**
 * Z is minimal polynomial t(x) = (x - 1)(x - 2)...(x - n)
 * s.t. it is equal to zero at all points that correspond to our gates.
 * Any polynomial that is equal to zero at all of these points must be
 * a multiple of Z.
 * If a polynomial is a multiple of Z, then its evaluation at any of those
 * points will be zero. We will exploit this equivalence.
 * Generate t(x) = (x - 1)(x - 2)...(x - n)
 *
 * @param degree The degree of the polynomial
 */
export function tPoly(degree: number): Fr[] {
    let t = [Fr.one()]
    for (let i = 1n; i <= degree; i++) {
        t = mulPolys(t, [new Fr(i).neg(), Fr.one()])
    }
    return t
}

/**
 * Quadratic arithmetic programme
 */
export interface QAP {
    /** Number of public inputs in the witness vector */
    l: number
    /** Lagrange basis for R1CS "A" gates */
    u: FrMatrix
    /** Lagrange basis for R1CS "B" gates */
    v: FrMatrix
    /** Lagrange basis for R1CS "C" gates */
    w: FrMatrix
}

/**
 * Create a quadratic arithmetic programme from an R1CS matrix triple (A, B, C)
 */
export function fromR1csToQap(l: number, A: FrMatrix, B: FrMatrix, C: FrMatrix): QAP {
    const A_t = transpose(A)
    const B_t = transpose(B)
    const C_t = transpose(C)

    // These are the x values or "positions" where each element of the transposed
    // matrices are evaluated at.
    const xs = Array.from({ length: A_t[0].length }, (_, i) => new Fr(BigInt(i + 1)))

    // Each element is a polynomial evaluation at its respective index.
    const Upolys = A_t.map((row) => interpolate(xs, row))
    const Vpolys = B_t.map((row) => interpolate(xs, row))
    const Wpolys = C_t.map((row) => interpolate(xs, row))

    return {
        l,
        u: Upolys,
        v: Vpolys,
        w: Wpolys,
    }
}

/**
 * Compute solution polynomials given a witness vector and a QAP.
 *
 * @param witness witness vector
 * @param qap Quadratic arithmetic programme
 * @returns Solution polynomials
 */
export function computeSolutionPoly(witness: Fr[], { u, v, w }: QAP) {
    const au = u.reduce((acc, row, i) => addPolys(acc, mulPolys([witness[i]], row)), [])
    const av = v.reduce((acc, row, i) => addPolys(acc, mulPolys([witness[i]], row)), [])
    const aw = w.reduce((acc, row, i) => addPolys(acc, mulPolys([witness[i]], row)), [])
    // h(x)*t(x)
    const ht = subPolys(mulPolys(au, av), aw)

    for (let i = 1n; i <= u[0].length; i++) {
        if (!evalPoly(ht, new Fr(i)).equals(Fr.zero())) {
            throw new Error(`Solution polynomial does not evaluate to zero at ${i}`)
        }
    }

    return {
        au,
        av,
        ht,
    }
}

/**
 * Divide the "solution polynomial" by the minimal polynomial t(x), and verify
 * that the remainder is zero.
 * The output here is h(x).
 *
 * @param ht h(x)t(x) from the solution polynomial
 * @param t minimal polynomial t(x)
 */
export function computeDivisorPoly(ht: Fr[], t: Fr[]): Fr[] {
    const [h, remainder] = divPolys(ht, t)
    for (const x in remainder) {
        if (!remainder[x].equals(Fr.zero())) {
            throw new Error(`Remainder is not zero at ${x}`)
        }
    }
    return h
}
