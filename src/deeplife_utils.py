mapping = {'Aba': 'A', 'Ace': 'X', 'Acr': 'X', 'Ala': 'A', 'Aly': 'K', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cas': 'C',
           'Ccs': 'C', 'Cme': 'C', 'Csd': 'C', 'Cso': 'C', 'Csx': 'C', 'Cys': 'C', 'Dal': 'A', 'Dbb': 'T', 'Dbu': 'T',
           'Dha': 'S', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'Glz': 'G', 'His': 'H', 'Hse': 'S', 'Ile': 'I', 'Leu': 'L',
           'Llp': 'K', 'Lys': 'K', 'Men': 'N', 'Met': 'M', 'Mly': 'K', 'Mse': 'M', 'Nh2': 'X', 'Nle': 'L', 'Ocs': 'C',
           'Pca': 'E', 'Phe': 'F', 'Pro': 'P', 'Ptr': 'Y', 'Sep': 'S', 'Ser': 'S', 'Thr': 'T', 'Tih': 'A', 'Tpo': 'T',
           'Trp': 'W', 'Tyr': 'Y', 'Unk': 'X', 'Val': 'V', 'Ycm': 'C', 'Sec': 'U', 'Pyl': 'O', 'Mhs': 'H', 'Snm': 'S',
           'Mis': 'S', 'Seb': 'S', 'Hic': 'H', 'Fme': 'M', 'Asb': 'D', 'Sah': 'C', 'Smc': 'C', 'Tpq': 'Y', 'Onl': 'X',
           'Tox': 'W', '5x8': 'X', 'Ddz': 'A'}


def three_to_one(three_letter_code):
    if three_letter_code[0].upper() + three_letter_code[1:].lower() not in mapping:
        return 'X'
    return mapping[three_letter_code[0].upper() + three_letter_code[1:].lower()]
