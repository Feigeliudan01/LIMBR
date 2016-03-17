os.chdir('/Volumes/jtkjtk/Jen/Jen/HTSeq')


def read_names(fname):
    paths = []
    with open(fname,'rU') as f:
        reader=csv.reader(f,delimiter='\t')
        for path in reader:
            paths.append(*path)
    return paths

def read_files(files):
    time = re.search('set(\d-.*)hours', files[0]).group(1)
    data = pd.read_csv(files[0], sep = '\t', header=None)
    data.columns = ['ind',time]
    data = data.set_index('ind')
    for fname in files[1:]:
        time = re.search('set(\d-.*)hours', fname).group(1)
        temp = pd.read_csv(fname, sep = '\t', header=None)
        temp.columns = ['ind',time]
        temp = temp.set_index('ind')
        data = pd.concat([data, temp], axis=1)
    return data


paths = read_names('names.txt')

data = read_files(paths)

data
