from ppiref.utils.misc import get_partition


def test_get_partition():
    # Test data
    lst = list(map(str, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
    partition_count = 5
    partition_size = 1 / partition_count

    # Perform the partition
    partitions = []
    for i in range(partition_count):
        beg = i * partition_size
        end = (i + 1) * partition_size
        partition = get_partition(lst, beg, end)
        partitions.append(partition)

    # Check disjointness and completeness
    flattened_partitions = [item for sublist in partitions for item in sublist]
    flattened_unique = list(set(flattened_partitions))
    assert len(flattened_partitions) == len(flattened_unique), "Partitions are not disjoint"
    assert sorted(flattened_partitions) == sorted(lst), "Partitions are not complete"

    # Test trivial partitioning
    assert get_partition(lst, 0., 1.) == lst
