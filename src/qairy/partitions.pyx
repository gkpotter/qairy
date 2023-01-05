from sage.misc.all import cached_function


def partitions_with_min(n, k, min_val=0, depth=0):
    if k == depth + 1 and n - sum(()) >= min_val:
        return ((n - sum(()),),)
    return (
        item + (i,)
        for i in [x for x in range(min_val, n + 1)]
        for item in partitions_with_min(n - i, k, min_val, depth=depth + 1)
    )


def decreasing_partitions_with_min(n, k, min_val=0):
    if k < 1:
        return
    if k == 1:
        if n >= min_val:
            yield (n,)
        return
    for i in range(min_val, n + 1):
        for result in decreasing_partitions_with_min(n - i, k - 1, min_val=i):
            yield result + (i,)


@cached_function
def cached_decreasing_red_partitions_with_min(n, k, r, min_val=0):
    return list(decreasing_red_partitions_with_min(n, k, r, min_val=0))


def decreasing_red_partitions_with_min(n, k, r, min_val=0):
    if k < 1:
        return
    if k == 1:
        if n >= min_val and n % r != 0:
            yield (n,)
        return
    for i in [x for x in range(min_val, n + 1) if x % r != 0]:
        for result in decreasing_red_partitions_with_min(n - i, k - 1, r, min_val=i):
            yield result + (i,)


def increasing_partitions_with_min(n, k, min_val=0):
    if k < 1:
        return
    if k == 1:
        if n >= min_val:
            yield (n,)
        return
    for i in range(min_val, n + 1):
        for result in increasing_partitions_with_min(n - i, k - 1, min_val=i):
            yield (i,) + result


@cached_function
def cached_increasing_red_partitions_with_min(n, k, r, min_val=0):
    return list(increasing_red_partitions_with_min(n, k, r, min_val=0))


def increasing_red_partitions_with_min(n, k, r, min_val=0):
    if k < 1:
        return
    if k == 1:
        if n >= min_val and n % r != 0:
            yield (n,)
        return
    for i in [x for x in range(min_val, n + 1) if x % r != 0]:
        for result in increasing_red_partitions_with_min(n - i, k - 1, r, min_val=i):
            yield (i,) + result


def tuple_automorphism_number(t):
    a = 1
    seen = []
    for x in t:
        seen.append(x)
        a *= seen.count(x)
    return a


# Find all pairs I,J whose union is s
def get_pairs(S):
    if len(S) == 0:
        return [((), ())]

    elif len(S) == 1:
        return [(S, ()), ((), S)]

    else:
        unordered_pairs = subsets_k(S, 2)
        pairs = []

        for pair in unordered_pairs:
            pairs.append(pair)

            if pair[1] != pair[0]:
                pairs.append((pair[1], pair[0]))

        pairs.extend([(S, ()), ((), S)])
        return pairs


def subsets_k(collection, k):
    yield from partition_k(collection, k, k)


def partition_k(collection, m, k):
    if len(collection) == 1:
        yield (collection,)
        return

    first = collection[0]
    for smaller in partition_k(collection[1:], m - 1, k):
        if len(smaller) > k:
            continue

        if len(smaller) >= m:
            for n, subset in enumerate(smaller):
                yield smaller[:n] + ((first,) + subset,) + smaller[n + 1 :]

        if len(smaller) < k:
            yield ((first,),) + smaller
