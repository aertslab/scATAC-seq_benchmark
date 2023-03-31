class ItemNotCachedException(Exception):
    pass


class LRUCache:
    def __init__(self, n):
        self.n = n
        self.keys = list()
        self.key2element = dict()

    def __contains__(self, key):
        return key in self.key2element.keys()

    def __getitem__(self, key):
        if key in self:
            self.keys.remove(key)
            self.keys.insert(0, key)
            return self.key2element[key]
        else:
            raise ItemNotCachedException

    def __setitem__(self, key, value):
        if len(self.keys) >= self.n:
            last_key = self.keys.pop()
            del self.key2element[last_key]
        if key not in self.key2element.keys():
            self.keys.insert(0, key)
            self.key2element[key] = value
        else:
            self.keys.remove(key)
            self.keys.insert(0, key)
            self.key2element[key] = value

    def __len__(self):
        return len(self.keys)


def cache(n):
    def _cache_decorator(func):
        if n <= 0:
            def wrapper_func(*args):
                return func(*args)
        else:
            cache = LRUCache(n)

            def wrapper_func(*args):
                if args in cache:
                    return cache[args]
                else:
                    result = func(*args)
                    cache[args] = result
                    return result
        return wrapper_func

    return _cache_decorator
