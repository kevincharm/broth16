import { expect } from 'chai'
import { setup, prove, verify } from '../src/groth16'
import { FrMatrix, fromR1csToQap } from '../src/qap'
import { Fr } from '@kevincharm/blstoise'

function prepareMatrix(matrix: number[][]): FrMatrix {
    return matrix.map((row) => row.map((x) => new Fr(BigInt(x))))
}

// Vitalik's famous example from https://medium.com/@VitalikButerin/quadratic-arithmetic-programs-from-zero-to-hero-f6d558cea649
// We want to prove that we know the solution to a cubic equation:
//  x^3 + x + 5 = 35 (where the answer x = 3).
// We define the equation as follows:
//  y = x**3
//  output = y + x + 5
//
// We flatten the equation to quadratic ones:
// t1 = x * x
// y = t1 * x
// t2 = y + x
// output = t2 + 5
//
// Note that t1 and t2 are *intermediate variables* that were created in order to decompose the original
// equation into a sequence of (at most) quadratic equations.
//
// Now we convert the gates (i.e. the multiplications and additions) into R1CS.
// In R1CS, we have groups of 3 vectors (a,b,c), and the solution to an R1CS programme is a witness s vector
// that satisfies the equation s.a * s.b = s.c (. = dot product)
//
// Mapping:
//      [1, x, output, t1, y, t2]
// The solution vector will consist of assignments for all of the above variables in that order.
//
// (a,b,c) for the 1st gate (t1 = x * x):
// a = [0, 1, 0, 0, 0, 0]
// b = [0, 1, 0, 0, 0, 0]
// c = [0, 0, 0, 1, 0, 0]
//
// (a,b,c) for the 2nd gate (y = t1 * x):
// a = [0, 0, 0, 1, 0, 0]
// b = [0, 1, 0, 0, 0, 0]
// c = [0, 0, 0, 0, 1, 0]
//
// (a,b,c) for the 3rd gate (t2 = y + x):
// a = [0, 1, 0, 0, 1, 0] <-- "y + x"
// b = [1, 0, 0, 0, 0, 0] <-- "addition check" of s.a
// c = [0, 0, 0, 0, 0, 1]
//
// (a,b,c) for the 4th gate (output = t2 + 5):
// a = [5, 0, 0, 0, 0, 1]
// b = [1, 0, 0, 0, 0, 0]
// c = [0, 0, 1, 0, 0, 0]

const R1CS = {
    A: prepareMatrix([
        [0, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 1, 0],
        [5, 0, 0, 0, 0, 1],
    ]),
    B: prepareMatrix([
        [0, 1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0],
    ]),
    C: prepareMatrix([
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0],
    ]),
}
const witness = [3, 35, 9, 27, 30].map((x) => new Fr(BigInt(x)))

describe('groth16', () => {
    it('compute and verify proof', () => {
        const l = 1
        const qap = fromR1csToQap(l, R1CS.A, R1CS.B, R1CS.C)
        const trustedSetup = setup(qap)
        const { proof, vk } = prove(trustedSetup, qap, witness)

        const publicInputs = witness.slice(0, l)
        expect(verify(proof, publicInputs, vk)).to.equal(true)
    }).timeout(10000)
})
