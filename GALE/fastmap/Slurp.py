"""
    This file is part of GALE,
    Copyright Joe Krall, 2014.

    GALE is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GALE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with GALE.  If not, see <http://www.gnu.org/licenses/>.
"""

from Row import *
from lib import *
from Table import *
  
def slurpFile(csvFile):
  t=want=None
  for _,cells in rows(csvFile):
    if t:
       if len(cells) == want: 
         t.put(t.read(cells)) 
    else:
       want, t = len(cells),  Table(cells)
  return t
  
def slurp(list, names):
    #Vivek
    #print list # [[decision + objectives (represented as '?')]] ? means missing values, it is a list of lists
    #print names  # [<$>decision.names (<<)objectives.names]


    t=want=None

    #print "Length of the list is : ", len(list)

    for rows in list:
        if t:
            if len(rows) == want: #check if it has the same numbers of columns as the first row
                t.put(rows)
        else:
            want, t = len(rows), Table(names) # it seems Joe has skipped the first row
            t.put(rows)  # added by Vivek

    # print "Number of rows in the table is: ", len(t.rows)
    # print "This is how a row looks like: ", t.rows[0]
    return t


