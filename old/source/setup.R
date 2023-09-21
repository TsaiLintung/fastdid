require(data.table)
require(stringr)
require(fixest)
require(rlist)
require(kit)
require(collapse)
require(dreamerr)


setDTthreads(0)
options(kit.nThread = getDTthreads())
setFixest_nthreads(getDTthreads())
