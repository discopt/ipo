import IPO
cimport cppIPO
import sage.all

#vectors
from sage.modules.free_module_element import vector as sagevec
#fields (usage: QQ = RationalField() vector(QQ, [...]))
from sage.rings.rational_field import RationalField
from sage.rings.integer_ring import IntegerRing

import rpy2.robjects as robjects
####################################
#Convert IPO-Vector to Sage Vector

def IPOVector_to_SageVector(IPOVector vec)
    #determine Field
    QQ = RationalField()
    #create list from IPO data
    sagelist = []
    size = vec.size()
    #find max index of the entries
    maxindex = 0
    for i in range(0,size):
        index = vec.index(i)
        if(index > maxindex):
            maxindex = index

    #initialize sagelist with zeros
    for (i in range(0, maxindex)):
        sagelist[i] = 0

    for (p in range(0,size)):
        index = vec.index(p)
        value = vec.value(index)
        sagelist[index] = value

    sagevector = sagevec(QQ, sagelist)
