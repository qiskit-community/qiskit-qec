import logging

import numpy as np

logger = logging.getLogger(__name__)


def generate_swath_error(
    codeword_array: np.ndarray,
    width: int,
    offset: int = 0,
    probability_of_error=1,
    is_vertical=True,
):
    """Expects ALL codewords to be numpy array of shape: (L//2, L);
    can achieve iid noise by setting width to width of codeword"""

    if probability_of_error <= 0:
        logger.warning("Did you mean to set error probability to 0?")
    error_mask = np.zeros(codeword_array.shape, dtype=int)

    num_rows = len(codeword_array)
    num_col = len(codeword_array[0])
    if is_vertical:
        if probability_of_error == 1:
            error_mask[:, offset : width + offset] = 1
        else:
            for i in range(num_rows):
                for j in range(width):
                    error_mask[i][(j + offset) % num_col] = np.random.choice(
                        [0, 1], p=[1 - probability_of_error, probability_of_error]
                    )

    else:
        height = width
        if probability_of_error == 1:
            error_mask[offset : offset + height, :] = 1
        else:
            for i in range(width):
                for j in range(num_col):
                    error_mask[(i + offset) % num_rows][j] = np.random.choice(
                        [0, 1], p=[1 - probability_of_error, probability_of_error]
                    )

    return error_mask ^ codeword_array, error_mask
