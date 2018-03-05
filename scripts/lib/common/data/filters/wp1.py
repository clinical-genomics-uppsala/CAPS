def evaluate_tuple(columns,mapper,condition):
    """
    """
    if isinstance(condition, tuple):
        return condition[0](columns,mapper,condition[1],condition[2])
    else:
        return condition(columns,mapper)

def evalute_else(columns,mapper,condition):
    """
    """
    if isinstance(condition, tuple) and len(condition) == 3:
        if condition[0](columns,mapper):
            return evaluate_tuple(columns,mapper,condition[1])
        else:
            return evaluate_tuple(columns,mapper,condition[2])
    else:
        raise Exception("Expecting a tuple of length 3")

def and_condition(columns,mapper,condition1,condition2):
    """
    """
    return evaluate_tuple(columns,mapper,condition1) and \
           evaluate_tuple(columns,mapper,condition2)

def or_condition(columns,mapper,condition1,condition2):
    """
    """
    return evaluate_tuple(columns,mapper,condition1) or \
           evaluate_tuple(columns,mapper,condition2)

def gene_name(name):
    """
    """
    return lambda columns, mapper: name == columns[mapper['Gene']]

def base_in_ref_or_var(base):
    """
    """
    return lambda columns, mapper: base in columns[mapper['Reference_base']] or base in columns[mapper['Variant_base']]

def vaf_above_or_equal(vaf):
    """
    """
    return lambda columns, mapper: float(columns[mapper['Variant_allele_ratio']]) >= vaf

def insertion_var(length):
    """
    """
    return lambda columns, mapper: int(columns[mapper['End']]) - int(columns[mapper['Start']]) == 0 and len(columns[mapper['Variant_base']]) == length

if __name__ == "__main__":
    import doctest
    doctest.testmod()
