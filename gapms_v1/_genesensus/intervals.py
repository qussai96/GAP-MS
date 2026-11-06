"""Small deterministic interval primitives and a dynamic coding-overlap index."""

from collections import defaultdict


def merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(b, out[-1][1]))
        else:
            out.append((a, b))
    return tuple(out)


def overlap(a, b):
    i = j = 0
    while i < len(a) and j < len(b):
        if max(a[i][0], b[j][0]) < min(a[i][1], b[j][1]):
            return True
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return False


def intersection_length(a, b):
    i = j = 0
    n = 0
    while i < len(a) and j < len(b):
        n += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def subtract(a, b):
    out = []
    b = merge(b)
    for left, right in merge(a):
        cursor = left
        for x, y in b:
            if y <= cursor:
                continue
            if x >= right:
                break
            if x > cursor:
                out.append((cursor, x))
            cursor = max(cursor, y)
            if cursor >= right:
                break
        if cursor < right:
            out.append((cursor, right))
    return tuple(out)


def disputed(a, b, window=9):
    regions = list(subtract(a.cds, b.cds)) + list(subtract(b.cds, a.cds))
    for x, y in set(a.introns) ^ set(b.introns):
        regions.extend(
            ((max(0, x - window), x + window), (max(0, y - window), y + window))
        )
    return merge(regions)


class CodingIndex:
    def __init__(self, bin_size=10000):
        self.bin_size = bin_size
        self.bins = defaultdict(set)
        self.items = {}

    def keys(self, model):
        return {
            (model.contig, model.strand, i)
            for a, b in model.cds
            for i in range(a // self.bin_size, (b - 1) // self.bin_size + 1)
        }

    def add(self, c):
        if c.uid in self.items:
            return
        self.items[c.uid] = c
        for k in self.keys(c.model):
            self.bins[k].add(c.uid)

    def remove(self, c):
        del self.items[c.uid]
        for k in self.keys(c.model):
            self.bins[k].discard(c.uid)

    def conflicts(self, c):
        ids = set()
        for k in self.keys(c.model):
            ids.update(self.bins[k])
        return [
            self.items[i]
            for i in sorted(ids)
            if i != c.uid and overlap(c.model.cds, self.items[i].model.cds)
        ]
