import re

try:
    from hgvslib.pHGVS import pHGVS
    from hgvs.easy import parse
except ModuleNotFoundError:
    pHGVS = None
    parse = None


AA_SUB_RE = re.compile(r"^[A-Z*]\d+(?:[A-Z*=]|Ter)$", re.IGNORECASE)
AA_FS_RE = re.compile(r"^[A-Z]\d+(?:[A-Z])?fs(?:\*(?:\d+)?|Ter\d*|X\d*)?$", re.IGNORECASE)
AA_DEL_DUP_RE = re.compile(r"^[A-Z]\d+(?:_[A-Z]\d+)?(?:del|dup)$", re.IGNORECASE)
AA_DELINS_RE = re.compile(r"^[A-Z](\d+)(?:_[A-Z](\d+))?delins([A-Z*]+)$", re.IGNORECASE)
AA_RANGE_INS_RE = re.compile(r"^[A-Z]\d+_[A-Z]\d+ins[A-Z*]+$", re.IGNORECASE)


def _is_aa_mut_light(alteration):
    text = str(alteration or "")
    return bool(
        AA_SUB_RE.match(text)
        or AA_FS_RE.match(text)
        or AA_DEL_DUP_RE.match(text)
        or AA_DELINS_RE.match(text)
        or AA_RANGE_INS_RE.match(text)
    )

def is_aa_mut_(alteration):
    if pHGVS is None:
        return is_aa_mut(alteration)
    type = pHGVS("p." + alteration).type
    if type != "?":
        return True
    else:
        return None


def is_aa_mut(alteration, transcript="NM_000000.00"):
    if _is_aa_mut_light(alteration):
        return True
    if parse is None:
        return False
    try:
        parse_aa_mut(alteration, transcript)
        return True
    except:
        return False


def parse_aa_mut(alteration, transcript= "NM_000000.00"):
    if parse is None:
        raise ImportError("hgvs is not installed")
    if alteration.startswith('p.'):
        variant_string = f"{transcript}:{alteration}"
    else:
        variant_string = f"{transcript}:p.{alteration}"
    var = parse(variant_string)
    return var


def get_aa_mut_attr(alteration):
    var = parse_aa_mut(alteration)
    var_dict = {
        # "variant": alteration,
        "start.aa": getattr(var.posedit.pos.start, 'aa', None),
        "start.pos": getattr(var.posedit.pos.start, 'pos', None),
        "end.aa": getattr(var.posedit.pos.end, 'aa', None),
        "end.pos": getattr(var.posedit.pos.end, 'pos', None),
        "edit.alt": getattr(var.posedit.edit, 'alt', None),
        "edit.length": getattr(var.posedit.edit, 'length', None),
        "edit.ref": getattr(var.posedit.edit, 'ref', None),
        "edit.type": getattr(var.posedit.edit, 'type', None),
    }
    return var_dict

def get_aa_mut_type(alteration):
    text = str(alteration or "")
    if AA_SUB_RE.match(text):
        if text.endswith("="):
            return "synonymous"
        if text.startswith("*"):
            return "stop lost"
        if "*" in text or text.lower().endswith("ter"):
            return "nonsense"
        return "missense"
    if AA_FS_RE.match(text):
        return "frameshift"
    if AA_DEL_DUP_RE.match(text):
        return "duplication" if text.lower().endswith("dup") else "in-frame deletion"
    if AA_RANGE_INS_RE.match(text):
        return "in-frame insertion"
    match = AA_DELINS_RE.match(text)
    if match:
        start = int(match.group(1))
        end = int(match.group(2) or start)
        return "in-frame insertion" if len(match.group(3)) > end - start + 1 else "in-frame deletion"

    if parse is None:
        return None

    try:
        attr = get_aa_mut_attr(alteration)
    except Exception:
        return None
    edit_type = attr["edit.type"]
    alt = attr["edit.alt"]
    ref = attr["edit.ref"]
    ref_s = str(ref) if ref is not None else ""
    alt_s = str(alt) if alt is not None else ""

    var_type = None
    stop_terms = {"ter", "*", "stop"}

    if edit_type == "sub":
        if alt_s in stop_terms:
            var_type = "nonsense"
        elif ref_s in stop_terms and alt_s not in stop_terms:
            var_type = "stop lost"
        elif ref_s == alt_s and ref_s != "":
            var_type = "synonymous"
        else:
            var_type = "missense"

    if edit_type == "fs":
        var_type = "frameshift"
    if edit_type == "dup":
        var_type = "duplication"
    if edit_type == "del":
        var_type = "in-frame deletion"
    if edit_type == "ins":
        var_type = "in-frame insertion"

    if edit_type == "delins":
        del_len = len(ref) if ref else 1
        ins_len = len(alt) if alt else 0
        if ins_len > del_len:
            var_type = "in-frame insertion"
        else:
            var_type = "in-frame deletion"

    if edit_type == "identity":
        var_type = "synonymous"

    # if var_type is None:
    #     print(alteration)

    return var_type


def is_cds_mut(alteration, transcript="NM_000000.00"):
    if parse is None:
        return bool(re.match(r"^c\.", alteration, re.IGNORECASE))
    try:
        parse_cds_mut(alteration, transcript)
        return True
    except:
        return False

def parse_cds_mut(alteration, transcript= "NM_000000.00"):
    if parse is None:
        raise ImportError("hgvs is not installed")
    # if alteration.startswith('c.'):
    #     variant_string = f"{transcript}:{alteration}"
    # else:
    #     variant_string = f"{transcript}:c.{alteration}"
    variant_string = f"{transcript}:{alteration}"
    var = parse(variant_string)
    return var

def get_cds_mut_attr(alteration):
    var = parse_cds_mut(alteration)
    var_dict = {
        "start.base": getattr(var.posedit.pos.start, 'base', None),
        "start.offset": getattr(var.posedit.pos.start, 'offset', None),
        "end.base": getattr(var.posedit.pos.end, 'base', None),
        "end.offset": getattr(var.posedit.pos.end, 'offset', None),
        "edit.alt": getattr(var.posedit.edit, 'alt', None),
        "edit.ref": getattr(var.posedit.edit, 'ref', None),
        "edit.type": getattr(var.posedit.edit, 'type', None),
    }
    return var_dict

def get_cds_mut_type(alteration):
    if parse is None:
        return "splice" if re.search(r"[+-]\d+", alteration) else "mutation"
    attr = get_cds_mut_attr(alteration)
    var_type = attr["edit.type"]
    off_s = attr["start.offset"]
    off_e = attr["end.offset"]
    off = off_s if off_s not in (None, 0) else off_e
    if off not in (None, 0):
        return "splice"
    return var_type






if __name__ == '__main__':
    # print(gene_norm("MLF1"))

    # print(type_norm('intragenic'))
    #
    # var_c = parse_cds_mut("c.790+1G>A")  #
    # var_c = parse_cds_mut("c.464-1_464-2delinsAT")
    # cc = get_aa_mut_attr("c.790+1G>A")
    # print(is_cds_mut('c.790+1G>A'))
    # print(is_aa_mut("D335Gfs"))
    # # var_g = parse("NM_000017.11:p.deletion")
    #
    # entity = 'deletion fff'
    # res = type_norm(entity)
    #
    #
    # res = get_aa_mut_type('E55=')
    res = get_aa_mut_type('R3128*')
    print(res)
    print('ok')
