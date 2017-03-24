from collections import defaultdict

def read_submission(input_fn):
    submission_dict = defaultdict(float)
    with open(input_fn, 'r') as input_file:
        input_file.readline()
        for line in input_file:
            key, value = line.strip().split(':')
            submission_dict[key] = float(value)
    return submission_dict


def cosine_similarity(d1, d2):
    d1_mag = sum([v**2 for v in d1.values()])
    d2_mag = sum([v**2 for v in d2.values()])
    dot_product = sum([d1[k]*d2[k] for k in d1.keys()])

    return dot_product/(d1_mag**.5 * d2_mag**.5)


def main(ans_fn, submission_fn):
    ans_dict = read_submission(ans_fn)
    submission_dict = read_submission(submission_fn)
    return 100*cosine_similarity(ans_dict, submission_dict)

if __name__ == "__main__":
    test_folder = 'hw5_W_0'
    ans_fn = './{}/{}_ans.txt'.format(test_folder, test_folder)
    submission_fn = './{}/{}_ans.txt'.format(test_folder, test_folder)
    main(ans_fn, submission_fn)
