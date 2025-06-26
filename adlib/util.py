import os
import re


def get_metal_from_path(path):
    m1 = re.search('/metal/([a-zA-Z]+)/', path)
    if m1:
        return m1[1]
    