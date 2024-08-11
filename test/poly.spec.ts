import { expect } from 'chai'
import { addPolys, divPolys, interpolate, mulPolys, subPolys, toPoly, transpose } from '../src/poly'
import { Fr } from '@kevincharm/blstoise'

function normalise(poly: (Fr | number | bigint)[]): string[] {
    return poly.map((x) => {
        if (x instanceof Fr) {
            return String(x.value)
        } else {
            return String(x)
        }
    })
}

describe('poly', () => {
    it('add', () => {
        const lhs = toPoly([1, 2, 3])
        const rhs = toPoly([3, 2, 1])
        expect(normalise(addPolys(lhs, rhs))).to.deep.eq(normalise([4, 4, 4]))
    })

    it('sub', () => {
        const lhs = toPoly([0, 1, 2])
        const rhs = toPoly([1, 1, 1])
        expect(normalise(subPolys(lhs, rhs))).to.deep.eq(
            normalise([Fr.zero().sub(new Fr(1n)).value, 0, 1]),
        )
    })

    it('mul', () => {
        const lhs = toPoly([1, 2, 3])
        const rhs = toPoly([7, 8, 9])
        expect(normalise(mulPolys(lhs, rhs))).to.deep.eq(
            normalise([1 * 7, 1 * 8 + 2 * 7, 1 * 9 + 2 * 8 + 3 * 7, 2 * 9 + 3 * 8, 3 * 9]),
        )
    })

    it('div', () => {
        const lhs = toPoly([1, 4])
        const rhs = toPoly([1, 2])
        const [quotient, remainder] = divPolys(lhs, rhs)
        expect(normalise(quotient)).to.deep.eq(normalise([2]))
        expect(normalise(remainder)).to.deep.eq(normalise([Fr.one().sub(new Fr(2n))]))
    })

    it('lagrange interpolation', () => {
        const l = interpolate(toPoly([1, 2, 3]), toPoly([4, 5, 6]))
        expect(normalise(l)).to.deep.eq(normalise([3, 1, 0]))
    })

    it('matrix transpose', () => {
        const transposed = transpose([
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ])
        expect(transposed).to.deep.eq([
            [1, 4, 7],
            [2, 5, 8],
            [3, 6, 9],
        ])
    })
})
