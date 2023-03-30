import pandas as pd
import numpy as np



# def flatten_columns(df, cols):
#     """Flattens multiple columns in a data frame, cannot specify all columns!"""
#     flattened_cols = {}
#     for col in cols:
#         flattened_cols[col] = pd.DataFrame([(index, value) for (index, values) in df[col].iteritems() for value in values.split(';')],
#                                            columns=['index', col]).set_index('index')
#     flattened_df = df.drop(cols, axis=1)
#     for col in cols:
#         flattened_df = flattened_df.join(flattened_cols[col])
#     return flattened_df

# df = pd.read_csv('escape_calculator_data.csv')

# encoding = (
#     df
#     [["eliciting_virus", "known_to_neutralize", "neg_log_IC50", "condition"]]
#     .drop_duplicates()
#     .reset_index(drop=True)
#     .assign(encoding=lambda x: x.index)
# )

# dms_data_encoded = (
#     df
#     .merge(encoding)
#     [["encoding", "site", "escape"]]
# )
# assert len(dms_data_encoded) == len(dms_data_encoded.drop_duplicates())

# encoding = encoding.drop(columns="condition")


# #
# #
# # CREATE NEW DF
# #
# #

# mut_rows = []

# knowns = ['BA.1', 'BA.2', 'BA.2.75', 'BA.5', 'BQ.1.1', 'D614G', 'XBB', 'any']
# viruses = ['SARS-CoV-1 then SARS-CoV-2', 'SARS-CoV-2', 'pre-Omicron SARS-CoV-2', 'pre-Omicron SARS-CoV-2 then BA.1', 'pre-Omicron SARS-CoV-2 then BA.2', 'pre-Omicron SARS-CoV-2 then BA.5']
# sites = [331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531]
# columns = None

# for param_known in knowns:
#     for virus in viruses:
#         for site in sites:
#             print(param_known, virus, site)
#             # Lookup
#             nencoding = encoding[['known_to_neutralize', 'neg_log_IC50', 'encoding']]
#             nencoding = pd.merge(dms_data_encoded, nencoding, on='encoding')

#             # Flatten
#             new_rows = []
#             new_index = []
#             new_cols = nencoding.columns
#             for index, row in nencoding.iterrows():
#                 neg_log_IC50 = row['neg_log_IC50'].split(';')
#                 known_to_neutralize = row['known_to_neutralize'].split(';')
#                 for n, k in zip(neg_log_IC50, known_to_neutralize):
#                     new_row = list(row[:3]) + [k, n]
#                     new_rows.append(new_row)
#                     new_index.append(index)
#             nencoding = pd.DataFrame(new_rows, columns=new_cols)

#             # Filter
#             nencoding = nencoding[nencoding['known_to_neutralize'] == param_known]

#             # Group by
#             nencoding = nencoding.groupby(['encoding', 'site', 'neg_log_IC50']).agg({'escape': 'mean'})
#             nencoding = nencoding.reset_index()

#             # Lookup
#             nnencoding = encoding[['eliciting_virus', 'encoding']]
#             nencoding = pd.merge(nencoding, nnencoding, on='encoding')

#             # Flatten
#             nencoding = flatten_columns(nencoding, ['eliciting_virus'])

#             # Filter
#             nencoding = nencoding[nencoding['eliciting_virus'] == virus]

#             # Group by
#             nencoding = nencoding.groupby(['encoding', 'site', 'neg_log_IC50']).agg({'escape': 'mean'})
#             nencoding = nencoding.reset_index()

#             # Group by apart
#             condition_escape_max = nencoding.groupby(['encoding']).agg({'escape': 'max'})
#             condition_escape_max = condition_escape_max.reset_index()
#             nencoding = pd.merge(nencoding, condition_escape_max, on='encoding')
#             nencoding['escape'] = nencoding['escape_x']
#             nencoding['condition_escape_max'] = nencoding['escape_y']
#             nencoding = nencoding.drop(columns=['escape_x', 'escape_y'])

#             # Transform
#             nencoding['site_binding_retained'] = 1
#             nencoding['site_binding_retained'] = nencoding['site_binding_retained'].where(nencoding['site'] != site, (nencoding['condition_escape_max'] - nencoding['escape']) / nencoding['condition_escape_max'])

#             # Transform
#             nencoding['encoding_weight'] = nencoding['neg_log_IC50']

#             # Group by apart
#             binding_retained = nencoding.groupby(['encoding', 'encoding_weight']).agg({'site_binding_retained': 'prod'})
#             binding_retained = binding_retained.reset_index()
#             nencoding = pd.merge(nencoding, binding_retained, on=['encoding','encoding_weight'], suffixes=['', '_new'])
#             nencoding['binding_retained'] = nencoding['site_binding_retained_new']
#             nencoding = nencoding.drop(columns=['site_binding_retained_new'])

#             # Transform
#             nencoding['encoding_weight'] = np.array(nencoding['encoding_weight'] , dtype='float32')
#             nencoding['escape_weighted'] = nencoding['encoding_weight'] * nencoding['escape']

#             # Transform
#             nencoding['escape_after_mut'] = np.power(nencoding['binding_retained'], 2) * nencoding['escape_weighted']
#             nencoding['binding_retained_exp'] = nencoding['encoding_weight'] * np.power(nencoding['binding_retained'], 2)

#             # Group by apart
#             n_conditions = float(nencoding.agg({'encoding': pd.Series.nunique}))

#             # Group by apart
#             new_trans = nencoding.groupby(['site']).agg({'escape_after_mut': 'sum', 'escape_weighted': 'sum'})
#             new_trans = new_trans.reset_index()
#             new_trans2 = nencoding.agg({'encoding_weight': 'sum', 'binding_retained_exp': 'sum'})
#             nencoding = pd.merge(nencoding, new_trans, on=['site'], suffixes=['', '_new'])
#             nencoding['sum_mutated'] = nencoding['escape_after_mut_new']
#             nencoding['sum_unmutated'] = nencoding['escape_weighted_new']
#             nencoding = nencoding.drop(columns=['escape_after_mut_new', 'escape_weighted_new'])

#             # Transform
#             nencoding['mutated'] = nencoding['sum_mutated'] / n_conditions
#             nencoding['unmutated'] = nencoding['sum_unmutated'] / n_conditions
#             nencoding['bound'] = new_trans2['binding_retained_exp'] / new_trans2['encoding_weight']
#             nencoding['escaped'] = 1 - nencoding['bound']

#             # Filter
#             nencoding = nencoding[nencoding['site'] == site]

#             # Clean
#             nencoding['count'] = nencoding.shape[0]
#             nencoding = nencoding[['mutated', 'unmutated', 'bound', 'escaped', 'count']]
#             nencoding = nencoding.drop_duplicates()
#             nencoding['param_known'] = param_known
#             nencoding['virus'] = virus
#             nencoding['site'] = site
            
#             if len(nencoding.values) != 0:
#                 mut_rows.append(nencoding.values[0])
            
#             if columns is None:
#                 columns = list(nencoding.columns)

# mut_df = pd.DataFrame(mut_rows, columns=columns)
# mut_df.to_csv('antigen.csv')


mut_df = pd.read_csv('escape_calculator.csv', index_col=0)
mut_df = mut_df[mut_df['virus'] == 'SARS-CoV-2']
mut_df = mut_df[mut_df['param_known'] == 'BA.1']
mut_df['synonymous'] = mut_df['escaped'] <= 0.15

print('Growth rate:', sum(mut_df['count'])/29904)
print(sum(mut_df['count']))
print(sum(mut_df[mut_df['synonymous'] == True]['count']))
print(sum(mut_df[mut_df['synonymous'] == False]['count']))

subset = mut_df[mut_df['synonymous'] == False]['count']

uniq = [1, sum(subset)]
n_mutations = 3

for i in range(2, 1 + n_mutations):
    cur = uniq[i-1] * uniq[1] 
    for k in range(i-1):
        cur -= (-1)**k * np.sum(np.power(subset, 2+k)) * uniq[i-2-k]
    cur //= i
    uniq.append(cur)

print(sum(mut_df['count'])**3)
print(sum(mut_df['count'])**3 * 25861)
print(uniq)
print(uniq[-1] * 25861)
print(sum(mut_df['count'])**3 * 25861 - uniq[-1] * 25861)

print(mut_df.head())