"""Genomic interval helpers."""

def merge_intervals(intervals):
    """
    intervals: [(start, end), ...], 1-based inclusive
    """
    if not intervals:
        return []

    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [list(intervals[0])]

    for s, e in intervals[1:]:
        if s <= merged[-1][1] + 1:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])

    return [(a, b) for a, b in merged]

def point_in_intervals(pos_1based, intervals_1based):
    for s, e in intervals_1based:
        if s <= pos_1based <= e:
            return True
    return False

def intervals_overlap(a_s, a_e, b_s, b_e):
    return not (a_e < b_s or b_e < a_s)

def any_overlap(a_s, a_e, intervals_1based):
    for s, e in intervals_1based:
        if intervals_overlap(a_s, a_e, s, e):
            return True
    return False