import string
import random


def random_id(length: int = 20) -> str:
    return ''.join(random.choices(string.digits, k=length))


def list_to_chunks(lst: list, n: int):
    """Cut list `lst` into `n` chunks.
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def get_partition(lst: list, beg: float, end: float):
    """Get [`beg`, `end`) partition of the list `lst`.
    """
    if not 0. <= beg <= end <= 1.:
        raise ValueError(
            "Invalid range. beg must be less than or equal to end,"
            "and both must be between 0 and 1."
        )

    n = len(lst)
    start_index = int(n * beg)
    end_index = int(n * end)

    return lst[start_index:end_index]
