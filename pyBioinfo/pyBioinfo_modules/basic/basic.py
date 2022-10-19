from functools import reduce


def makeList(ele1, ele2):
    if not isinstance(ele1, list):
        ele1 = [ele1]
    if not isinstance(ele2, list):
        ele2 = [ele2]
    return ele1 + ele2


def flattenList(lst):
    # print(flattenList([[[2, 3], [42, 'lskdj', 3]], [23, 8888, 88]]))
    if all(not isinstance(ele, list) for ele in lst):
        return lst
    return flattenList(reduce(makeList, lst, []))
