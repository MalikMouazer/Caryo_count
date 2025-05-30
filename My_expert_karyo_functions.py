import re
from collections import Counter
import pandas as pd

# Extraction des numéros de chromosome dans une anomalie ISCN
def get_chromosomes(anom):
    """
    Retourne l'ensemble des chromosomes impliqués dans l'anomalie,
    basé sur les notations avant chaque parenthèse de type der, del, dup, ins, t, i, ider, idic.
    """
    nums = set()
    for m in re.finditer(r'(?:der|del|dup|ins|t|i|ider|idic|r)\(([0-9;]+)', anom):
        for num in m.group(1).split(';'):
            nums.add(num)
    return nums

# Parsing de la formule karyotypique
def parse_caryotype(chaine_iscn):
    """
    Parse une chaîne ISCN avec clones séparés par '/'.
    Renvoie la liste plate des anomalies et un dict {anom: [clones]}.
    Gère aussi les anomalies de ploidie (≠46).
    """
    # Remove all whitespace for robust parsing
    chaine_iscn = re.sub(r"\s+", "", chaine_iscn)
    
    anomalies = []
    clone_map = {}
    clones = [re.sub(r"\[.*?\]", "", c) for c in chaine_iscn.split('/')]
    for idx, clone in enumerate(clones, start=1):
        parts = [p.strip().strip('.') for p in clone.split(',') if p.strip()]
        # Détection de la ploidie
        try:
            total = int(re.sub(r"\D", "", parts[0]))
            if total != 46:
                if total == 92:
                    pl = 'Tetraploidy'
                elif total == 69:
                    pl = 'Triploidy'

                anomalies.append(pl)
                clone_map.setdefault(pl, []).append(f"clone{idx}")
        except Exception:
            pass
        # Extraction des anomalies structurelles
        for an in parts[2:]:
            anomalies.append(an)
            clone_map.setdefault(an, []).append(f"clone{idx}")
    return anomalies, clone_map

# Détection des anomalies unichromosomiques déséquilibrées de poids 2
def is_single_chr_deseq(anom, count):
    """
    Détecte les anomalies unichromosomiques déséquilibrées qui valent 2 points:
    - Tetrasomie/triplication/quadruplication
    - Chromosome isodérivé
    """
    # Tetrasomie/triplication/quadruplication
    if anom.startswith('+') and count > 1:
        return True
    if anom.startswith('trp'):
        return True
    # Chromosome isodérivé ou isodicentrique
    if anom.startswith('ider'):
        return True
    return False

# Détection des anomalies équilibrées
def is_balanced_translocation(anom):
    """
    Détecte les translocations équilibrées:
    t(NUM;NUM[;...])(p;q) sans der,+,-
    """
    pattern = r'^t\(\d+(?:;\d+)+\)\(.+\)$'
    return bool(re.match(pattern, anom)) and 'der' not in anom and '+' not in anom and '-' not in anom

def is_unbalanced_translocation(anom):
    """
    Détecte les translocations déséquilibrées:
    - chromosome dérivé (der(...)) contenant un t(...) ou
    - tout t(...) non pure
    """
    # Cas d'un chromosome dérivé comportant une translocation
    if 'der' in anom and 't(' in anom:
        return True
    # Cas d'un t(...) quelconque non pur (équilibré)
    if 't(' in anom and not is_balanced_translocation(anom):
        return True
    return False

def is_balanced_insertion(anom):
    """
    Détecte les insertions équilibrées:
    ins(NUM;NUM[;...])(p;q1q2) sans der,+,-
    """
    pattern = r'^ins\(\d+(?:;\d+)+\)\(.+\)$'
    return bool(re.match(pattern, anom)) and 'der' not in anom and '+' not in anom and '-' not in anom

# Détection des anomalies multichromosomiques déséquilibrées pour 2 points
def is_complex_multichr_deseq(anom):
    """
    Détecte les anomalies multichromosomiques déséquilibrées (≥2 chromosomes) pour 2 points.
    Renvoie False si un seul chromosome impliqué.
    """
    chroms = get_chromosomes(anom)
    # si un seul chromosome impliqué -> pas multi-chromosomique déséquilibrée
    if len(chroms) <= 1:
        return False
    # dérivé ou anneau -> complexe multi-chromosomique
    if anom.startswith('der') or anom.startswith('r('):
        return True
    # insertion non pure -> complexe
    if 'ins(' in anom and not is_balanced_insertion(anom):
        return True
    # translocation non pure -> complexe
    if 't(' in anom and not is_balanced_translocation(anom):
        return True
    return False

# Typage pour affichage
def type_anomalie(anom):
    """
    Détermine le type d'anomalie pour l'affichage.
    Retourne une chaîne décrivant le type d'anomalie.
    """
    if is_complex_multichr_deseq(anom):
        return 'Multichromosomique déséquilibrée'
    if is_balanced_translocation(anom):
        return 'Translocation équilibrée'
    if is_unbalanced_translocation(anom):
        return 'Translocation déséquilibrée'
    if is_balanced_insertion(anom):
        return 'Insertion équilibrée'
    if anom == '<2n>':
        return 'Ploidy'
    if '~' in anom:
        return 'Pléiade chromosomique'
    if anom == '+mar':
        return 'Chromosome marqueur'
    if 'dmin' in anom:
        return 'Double minutes'
    if anom.startswith('hsr'):
        return 'Homogeneously staining region'
    if anom.startswith('r('):
        return 'Anneau'
    if anom.startswith('der'):
        return 'Chromosome dérivé'
    if anom.startswith('ins'):
        return 'Insertion'
    if anom.startswith('t('):
        return 'Translocation'
    if anom.startswith('+'):
        return 'Gain chr' + re.sub(r"\D", "", anom)
    if anom.startswith('-'):
        return 'Perte chr' + re.sub(r"\D", "", anom)
    if anom.startswith('dup'):
        return 'Duplication'
    if anom.startswith('del'):
        return 'Délétion'
    if anom.startswith('trp'):
        return 'Triplication/Quadruplication'
    if anom.startswith('idic'):
        return 'Isodicentric chromosome'
    if anom.startswith('ider'):
        return 'Isoderivative chromosome'
    if anom.startswith('i(') or 'iso' in anom:
        return 'Isochromosome'
    return 'Autre'

# Calcul des scores
def calcul_scores(anomalies, clone_map):
    """
    Calcule les scores selon deux méthodes:
    
    Jondreville 2020 : 1 point par anomalie (toujours)
    ISCN 2024     :
      - 0 pt pour anomalies implicites (dérivés implicites & constitutionnelles)
      - 2 pts pour déséquilibres unichr/multichr ou translocations déséquilibrées
      - 1 pt pour anomalies standard
    """
    counts = Counter(anomalies)

    # 1) Collecte des dérivés implicites pour chaque dérivé suivi d'un t(A;B)
    implicit = set()
    t_events = {}
    for an in counts:
        m = re.match(r"der\((\d+)\).*t\((\d+);(\d+)\)", an)
        if m:
            _, A, B = m.groups()
            key = tuple(sorted([A, B]))
            t_events.setdefault(key, []).append(an)
    for ders in t_events.values():
        # on considère « explicites » ceux qui ont add/del/dup
        explicits = [d for d in ders if re.search(r"add|del|dup", d)]
        if explicits:
            for d in ders:
                if d not in explicits:
                    implicit.add(d)

    # 2) Repérage des chr issus de dérivés multi-chr (der(X;Y))
    multi_der_chrs = set()
    for an in counts:
        if an.startswith("der"):
            m = re.match(r"^der\(([0-9;]+)\)", an)
            if m:
                chrs = m.group(1).split(";")
                if len(chrs) > 1:
                    multi_der_chrs.update(chrs)

    rows = []
    total_j = total_i = 0

    for anom, cnt in counts.items():
        score_j = 1  # Jondreville = 1 pour toutes

        # a) Constitutionnelles (+Nc) → ISCN = 0
        if re.match(r"^\+\d+c$", anom):
            score_i = 0
            explication = "Anomalie constitutionnelle (0 point)"

        # b) Dérivés implicites repérés plus haut → ISCN = 0
        elif anom in implicit:
            score_i = 0
            explication = "Dérivé implicite (0 point)"

        # c) Gains/pertes simples dans un der multi-chr → ISCN = 0
        elif anom.startswith(("+", "-")):
            num = re.sub(r"\D", "", anom)
            if num in multi_der_chrs:
                score_i = 0
                explication = "Gain/perte simple dans un dérivé multi-chromosomique (0 point)"
            else:
                if is_single_chr_deseq(anom, cnt):
                    score_i = 2
                    explication = "Déséquilibre unichromosomique (2 points)"
                elif is_complex_multichr_deseq(anom):
                    score_i = 2
                    explication = "Déséquilibre multichromosomique complexe (2 points)"
                elif is_unbalanced_translocation(anom):
                    score_i = 2
                    explication = "Translocation déséquilibrée (2 points)"
                else:
                    score_i = 1
                    explication = "Anomalie standard (1 point)"

        # d) Toutes les autres anomalies → scoring standard
        else:
            if is_single_chr_deseq(anom, cnt):
                score_i = 2
                explication = "Déséquilibre unichromosomique (2 points)"
            elif is_complex_multichr_deseq(anom):
                score_i = 2
                explication = "Déséquilibre multichromosomique complexe (2 points)"
            elif is_unbalanced_translocation(anom):
                score_i = 2
                explication = "Translocation déséquilibrée (2 points)"
            else:
                score_i = 1
                explication = "Anomalie standard (1 point)"

        total_j += score_j
        total_i += score_i

        rows.append({
            "Anomalie": anom,
            "Type": type_anomalie(anom),
            "Explication": explication,
            "Occurrences": cnt,
            "Clones": ", ".join(clone_map.get(anom, [])),
            "Score Jondreville 2020": score_j,
            "Score ISCN 2024": score_i,
        })

    # Ligne de totaux
    rows.append({
        "Anomalie": "TOTAL",
        "Type": "",
        "Explication": "",
        "Occurrences": "",
        "Clones": "",
        "Score Jondreville 2020": total_j,
        "Score ISCN 2024": total_i,
    })

    return pd.DataFrame(rows), total_i

# Fonction pour analyser une formule caryotypique
def analyser_formule(formule):
    """
    Analyse une formule caryotypique et retourne:
    - Le DataFrame des anomalies détectées
    - Le score total
    - Une erreur éventuelle
    """
    try:
        anomalies, clone_map = parse_caryotype(formule)
        df, total = calcul_scores(anomalies, clone_map)
        return df, total, None
    except Exception as e:
        return None, 0, f"Erreur lors de l'analyse de la formule: {str(e)}"
