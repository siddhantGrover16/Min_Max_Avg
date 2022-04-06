
from Bio import Phylo
import numpy as np
from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.BaseTree import Clade
from io import StringIO

from MinPD import GetMinPD
from MaxPD import GetMaxPD
from Avg import GetAvgPD


def getTreeFromFile(file):
    tree = Phylo.read(file,"newick")
    return tree

def getTreeFromString(string):
    tree = Phylo.read(StringIO(string), "newick")
    return tree




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print('PyCharm')


tree_string = "((((((((((A/swine/Colorado/A02245370/2019|H1N2|delta1a|1998|2019-12-27:5.1E-4,A/swine/Colorado/A02245372/2019|H1N2|TTTTPT|delta1a|1998|2019-12-31:0.00116)[&label=0.994]:0.0048,(A/swine/Colorado/A02245365/2019|H1N2|delta1a|1998|2019-12-19:0.00118,A/swine/Kansas/A02245678/2020|H1N2|delta1a|1998|2020-07-02:0.00297)[&label=0.961]:0.00238)[&label=0.811]:5.8E-4,A/swine/Colorado/A02245366/2019|H1N2|delta1a|1998|2019-12-19:0.00237)[&label=0.95]:0.00237,A/swine/Colorado/A01785480/2018|H1N2|delta1a|1998|2018-03-21:5.1E-4)[&label=1.0]:0.00834,A/swine/Colordao/A01785769/2018|H1N2|delta1a|1998|2018-11-28:0.01083)[&label=0.161]:5.1E-4,((A/swine/Colorado/A01785758/2018|H1N2|TTTTPT|delta1a|1998|2018-11-15:0.00116,A/swine/Colorado/A02245360/2019|H1N2|delta1a|1998|2019-12-18:0.0060)[&label=0.987]:0.0049,A/swine/Colorado/A01785307/2017|H1N2|TTTTPT|delta1a|1998|2017-10-03:0.00654)[&label=0.963]:0.00358)[&label=0.947]:0.0023,((A/swine/Iowa/A02157802/2018|H1N2|delta1a|2002|2018-04-12:0.01654,A/swine/Iowa/A02219785/2017|H1N2|TTPTPT|delta1a|2002|2017-08-24:0.00685)[&label=0.07]:5.2E-4,A/swine/Nebraska/A02172375/2018|H1N2|delta1a|2002|2018-07-12:0.01088)[&label=0.797]:0.0013)[&label=0.484]:5.1E-4,((((A/swine/Iowa/A02158150/2018|H1N2|delta1a|2002|2018-04-26:0.00298,A/swine/Iowa/A02218172/2017|H1N2|delta1a|2002|2017-06-26:5.1E-4)[&label=0.814]:5.9E-4,A/swine/Iowa/A02221492/2017|H1N2|delta1a|2002|2017-09-12:0.00119)[&label=0.661]:5.1E-4,A/swine/Iowa/A02431617/2019|H1N1|TTTPPT|delta1a|classicalSwine|2019-03-18:0.0072)[&label=0.993]:0.00475,A/swine/Iowa/A01668210/2016|H1N2|TTTTPT|delta1a|2002|2016-11-11:0.00119)[&label=0.753]:6.0E-4)[&label=0.918]:0.0028,((A/swine/Colorado/A02155470/2018|H1N2|TTTTPT|delta1a|2002|2018-03-05:0.00182,A/swine/Colorado/A02428130/2018|H1N2|TTTTPT|delta1a|2002|2018-12-27:0.00772)[&label=0.898]:0.00184,A/swine/Minnesota/A01785324/2017|H1N2|delta1a|2002|2017-10-20:0.00356)[&label=0.999]:0.0074)[&label=0.988]:0.00727,A/swine/Kansas/A01785288/2017|H1N2|delta1a|2002|2017-09-05:0.0139)[&label=0.989]"
tree_file = "h1n2.tre"
leaf_set_size = 15

MinPD = GetMinPD(getTreeFromString(tree_string) ,leaf_set_size)
MaxPD = GetMaxPD(getTreeFromString(tree_string) ,leaf_set_size)
AvgPD = GetAvgPD(getTreeFromString(tree_string) ,leaf_set_size)

print ("Min PD for given tree is ",MinPD )
print ("Max PD for given tree is ",MaxPD )
print ("Avg PD for given tree is ",AvgPD )



