"""Scientific intensity and full-grid regressions for gridrec."""

import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from tomopy.recon.algorithm import recon


class GridrecIntensityTestCase(unittest.TestCase):
    def test_blob_gain_power_of_two(self):
        self._check_blob_gain(256)

    def test_blob_gain_non_power_of_two(self):
        self._check_blob_gain(255)

    def _check_blob_gain(self, width):
        # Both detector widths use a 256-pixel internal FFT grid. Analytic
        # Gaussian projections avoid an interpolation-dependent projector.
        theta = np.linspace(0, np.pi, 450, endpoint=False, dtype=np.float32)
        detector = np.arange(width)
        sigma = 1.5
        positions = [(0, 0), (0, 0)]
        for fraction in (0.5, 0.85):
            offset = round(fraction * 128)
            for position in ((offset, 0), (-offset, 0),
                             (0, offset), (0, -offset)):
                # Every position is reconstructed in both native FFT channels.
                positions.extend((position, position))

        projections = []
        for x, y in positions:
            projected_center = x * np.cos(theta) + y * np.sin(theta) + width / 2
            projections.append(
                np.exp(-0.5 * ((detector[None, :] -
                               projected_center[:, None]) / sigma)**2) /
                (np.sqrt(2 * np.pi) * sigma))
        projections = np.stack(projections, axis=1).astype(np.float32)

        normalized = {}
        for algorithm in ('gridrec', 'fbp'):
            # One native call preserves pairs of real/imaginary slices.
            result = recon(projections, theta, center=width / 2,
                           algorithm=algorithm, filter_name='shepp', ncore=1)
            integrals = np.array([
                self._blob_integral(image, x, y)
                for image, (x, y) in zip(result, positions)
            ]).reshape(-1, 2)
            # Each FFT channel has its own central reference. Compare spatial
            # gain independently of the algorithms' absolute normalization.
            normalized[algorithm] = (integrals / integrals[0]).ravel()

        for index, (x, y) in enumerate(positions[2:], start=2):
            with self.subTest(width=width, x=x, y=y, channel=index % 2):
                assert_allclose(normalized['gridrec'][index],
                                normalized['fbp'][index], rtol=0.01, atol=0)

    @staticmethod
    def _blob_integral(image, x, y):
        half_width = 8
        # Recenter the integration window on the local peak to accommodate
        # the subpixel placement difference between gridrec and fbp.
        row = image.shape[0] // 2 - 1 - y
        column = image.shape[1] // 2 + x
        window = image[row - half_width:row + half_width + 1,
                       column - half_width:column + half_width + 1]
        peak_row, peak_column = np.unravel_index(
            np.argmax(window), window.shape)
        row += peak_row - half_width
        column += peak_column - half_width
        return image[row - half_width:row + half_width + 1,
                     column - half_width:column + half_width + 1].sum()


def _disk_projections(width):
    """Analytic line integrals through two off-axis absorbing disks."""
    theta = np.linspace(0, np.pi, 13, endpoint=False, dtype=np.float32)
    detector = np.arange(width, dtype=np.float32) - width / 2
    sino = np.zeros((len(theta), width), dtype=np.float32)
    for x, y, radius, attenuation in ((-1, .5, 4, .03), (2, -2, 1.8, .08)):
        distance = detector - (x * np.cos(theta) + y * np.sin(theta))[:, None]
        chord = 2 * np.sqrt(np.maximum(radius**2 - distance**2, 0))
        sino += attenuation * chord
    scale = 1 + np.arange(16, dtype=np.float32)[:, None, None] / 16
    return scale * sino, theta


def _reconstruct(data, theta, **kwargs):
    return recon(data, theta, algorithm='gridrec', sinogram_order=True,
                 ncore=1, filter_name='shepp', **kwargs)


class GridrecFullGridTestCase(unittest.TestCase):
    def test_outer_edges(self):
        # Samples from every outer edge, including all four corners, for the
        # analytic two-disk input with the corrected full-grid PSWF table.
        # Fixed values also detect a consistently wrong edge that would agree
        # between otherwise identical full-grid, crop, and batching calls.
        expected = {
            16: [[-0.014514812, -0.017745465, -0.011935343, -0.020524167],
                 [-0.090851687, -0.013721272, -0.006832186, -0.053797100],
                 [-0.014514812, -0.024034632, -0.011040418, -0.090851687],
                 [-0.020524167, -0.006192388, -0.011042548, -0.053797100]],
            15: [[-0.019327067, -0.012183571, -0.018951282, -0.018941890],
                 [-0.102695331, -0.029816074, -0.016673394, -0.033196878],
                 [-0.019327067, -0.039638739, -0.008431993, -0.102695331],
                 [-0.018941890, -0.000608895, -0.018098855, -0.033196878]],
        }
        for width in (16, 15):
            with self.subTest(width=width):
                data, theta = _disk_projections(width)
                full = _reconstruct(data, theta, num_gridx=16, num_gridy=16)
                edges = np.stack((full[:, 0, ::5], full[:, -1, ::5],
                                  full[:, ::5, 0], full[:, ::5, -1]), axis=1)
                scale = 1 + np.arange(16)[:, None, None] / 16
                assert_allclose(edges, scale * np.array(expected[width]),
                                rtol=1e-2, atol=1e-7)

    def test_full_grid_and_centered_crops(self):
        for width in (16, 15):
            data, theta = _disk_projections(width)
            full = _reconstruct(data, theta, num_gridx=16, num_gridy=16)
            self.assertTrue(np.isfinite(full).all())
            for nx, ny in ((16, 16), (15, 15), (14, 14), (16, 13),
                           (13, 16), (15, 12), (12, 15)):
                with self.subTest(width=width, grid=(nx, ny)):
                    crop = _reconstruct(data, theta, num_gridx=nx, num_gridy=ny)
                    self.assertEqual(crop.shape, (16, nx, ny))
                    # The output transposes the FFT grid and reverses its first
                    # axis, so odd crops start at ceil(dx/2), floor(dy/2).
                    row, col = (16 - nx + 1) // 2, (16 - ny) // 2
                    assert_array_equal(
                        crop, full[:, row:row + nx, col:col + ny])

    def test_full_grid_batching(self):
        for width in (16, 15):
            data, theta = _disk_projections(width)
            full = _reconstruct(data, theta, num_gridx=16, num_gridy=16)
            for batch in (5, 6, 16):
                with self.subTest(width=width, batch=batch):
                    pieces = [
                        _reconstruct(data[start:start + batch].copy(), theta,
                                     num_gridx=16, num_gridy=16)
                        for start in range(0, len(data), batch)
                    ]
                    result = np.concatenate(pieces)
                    # Odd batches change the native FFT partners. Use the
                    # existing gridrec fixture tolerance for float roundoff.
                    assert_allclose(result, full, rtol=1e-2, atol=1e-7)
