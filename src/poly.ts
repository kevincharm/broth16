import { Fr } from '@kevincharm/blstoise'

// Most of these functions come from Vitalik's QAP implementation:
// https://github.com/ethereum/research/blob/master/zksnark/qap_creator.py

/**
 * Helper to convert a list of coefficients to a list of polynomial
 * curve scalars coefficients.
 *
 * @param coeffs Coefficients
 * @returns Polynomial in F_r
 */
export function toPoly(coeffs: (number | bigint)[]): Fr[] {
    return coeffs.map((c) => new Fr(BigInt(c)))
}

/**
 * Create a zero poly
 *
 * @param degree Degree of polynomial
 * @returns Polynomial with all coefficients zero
 */
export function zeroPoly(degree: number): Fr[] {
    return Array.from({ length: degree }, () => Fr.zero())
}

/**
 * Add two polynomials
 *
 * @param lhs Addend
 * @param rhs Addend
 * @returns Sum
 */
export function addPolys(lhs: Fr[], rhs: Fr[]): Fr[] {
    const length = Math.max(lhs.length, rhs.length)
    const result: Fr[] = Array.from({ length }, () => Fr.zero())
    for (let i = 0; i < lhs.length; i++) {
        result[i] = result[i].add(lhs[i])
    }
    for (let i = 0; i < rhs.length; i++) {
        result[i] = result[i].add(rhs[i])
    }
    return result
}

/**
 * Subtract two polynomials
 *
 * @param lhs Minuend
 * @param rhs Subtrahend
 * @returns Difference
 */
export function subPolys(lhs: Fr[], rhs: Fr[]): Fr[] {
    const length = Math.max(lhs.length, rhs.length)
    const result: Fr[] = Array.from({ length }, (_, i) => (i < lhs.length ? lhs[i] : Fr.zero()))
    for (let i = 0; i < rhs.length; i++) {
        result[i] = result[i].sub(rhs[i])
    }
    return result
}

/**
 * Polynomial multiplication
 *
 * @param lhs Multiplicand
 * @param rhs Multiplier
 * @returns Multiplication
 */
export function mulPolys(lhs: Fr[], rhs: Fr[]): Fr[] {
    const length = lhs.length + rhs.length - 1
    const result: Fr[] = Array.from({ length }, () => Fr.zero())
    for (let i = 0; i < lhs.length; i++) {
        for (let j = 0; j < rhs.length; j++) {
            result[i + j] = result[i + j].add(lhs[i].mul(rhs[j]))
        }
    }
    return result
}

/**
 * Polynomial division
 *
 * @param lhs Dividend
 * @param rhs Divisor
 * @returns Quotient and remainder
 */
export function divPolys(lhs: Fr[], rhs: Fr[]): [quotient: Fr[], remainder: Fr[]] {
    const length = lhs.length - rhs.length + 1
    const result = Array.from({ length }, () => Fr.zero())
    let remainder = lhs // ref
    while (remainder.length >= rhs.length) {
        const leadingFactor = remainder.at(-1)!.mul(rhs.at(-1)!.inv())
        const pos = remainder.length - rhs.length
        result[pos] = leadingFactor
        const spent = [...zeroPoly(pos), leadingFactor]
        remainder = subPolys(remainder, mulPolys(rhs, spent)).slice(0, -1)
    }
    return [result, remainder]
}

/**
 * Evaluate a polynomial p(x)
 *
 * @param p Polynomial
 * @param x Position
 * @returns p(x)
 */
export function evalPoly(p: Fr[], x: Fr): Fr {
    return p.reduce((acc, coeff, i) => acc.add(coeff.mul(x.exp(BigInt(i)))), Fr.zero())
}

/**
 * Compute a polynomial p(x) s.t. p(x) = 0 foreach x in xs
 * Borrowed from https://github.com/GuildOfWeavers/galois/blob/master/lib/PrimeField.ts#L992
 *
 * @param xs
 * @returns
 */
function zpoly(xs: Fr[]): Fr[] {
    const result: Fr[] = Array.from({ length: xs.length + 1 }, () => Fr.zero())
    result[result.length - 1] = Fr.one()

    for (let i = 0, p = result.length - 2; i < xs.length; i++, p--) {
        result[p] = Fr.zero()
        for (let j = p; j < result.length - 1; j++) {
            result[j] = result[j].sub(result[j + 1].mul(xs[i]))
        }
    }
    return result
}

/**
 * Lagrange interpolation
 * Borrowed from https://github.com/GuildOfWeavers/galois/blob/f3d9cfbf2fe7857f3840bdba3406e2ba9ea548c7/lib/PrimeField.ts#L808
 *
 * @param xs
 * @param ys
 * @returns Interpolated polynomial fitting f(xs) = ys
 */
export function interpolate(xs: Fr[], ys: Fr[]): Fr[] {
    if (xs.length !== ys.length) {
        throw new Error('length mismatch')
    }

    const root = zpoly(xs)
    const divisor = [Fr.zero(), Fr.one()]
    const numerators: Fr[][] = Array.from({ length: xs.length }, () => [])
    for (let i = 0; i < xs.length; i++) {
        divisor[0] = xs[i].neg()
        ;[numerators[i]] = divPolys(root, divisor)
    }

    const denominators: Fr[] = Array.from({ length: xs.length }, () => Fr.zero())
    for (let i = 0; i < xs.length; i++) {
        denominators[i] = evalPoly(numerators[i], xs[i])
    }
    const invDenValues = denominators.map((d) => d.inv())

    const r = Array.from({ length: xs.length }, () => Fr.zero())
    for (let i = 0; i < xs.length; i++) {
        const ySlice = ys[i].mul(invDenValues[i])
        for (let j = 0; j < xs.length; j++) {
            if (!numerators[i][j].equals(Fr.zero()) && !ys[i].equals(Fr.zero())) {
                r[j] = r[j].add(numerators[i][j].mul(ySlice))
            }
        }
    }
    return r
}

/**
 * Transpose a matrix
 * Matrix must be correctly formed!
 *
 * e.g.
 * [[1, 2, 3],
 *  [4, 5, 6]]
 * =>
 * [[1, 4],
 *  [2, 5],
 *  [3, 6]]
 *
 * @param matrix
 * @returns Transposed matrix
 */
export function transpose<T>(matrix: T[][]) {
    // Assumes an m x n matrix is correctly formed
    const m = matrix.length // rows
    const n = matrix[0].length // cols
    const out: T[][] = Array.from({ length: n }, () => Array.from({ length: m }, () => null!))
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < m; j++) {
            out[i][j] = matrix[j][i]
        }
    }
    return out
}
