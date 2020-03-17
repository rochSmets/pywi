
import numpy as np


class Shape(object):

   """

   """

   def __init__(self,
                data,
                limit,
                chunk = None,
                stride = None,
                section = None):

      # ..... private members
      self._rank = data.ndim
      self._size = data.shape
      self._step = []

      # ..... public members
      self.axis = []
      self.data = np.empty(self._size)


      # ..... set the steps of the data's grid in each directions
      self._setStep(limit)

      # ..... set extract chunk from data with the given stride along each axes
      self._extractData(data, limit, chunk, stride)

      # ..... slice data along the given section
      self._sliceData(section)



   def _setStep(self,
                limit) :

      """
      this method set the grid step in each direction for the input data
      and the axis list from limit

      input : limit is a list of set [min, max] for each direction
              in physical unit
      """


      # set the steps in a list
      if limit.__len__() == self._rank :
         for i in range(self._rank) :
            step = (limit[i][1]-limit[i][0])/(self._size[i]-1.0)
            self._step.append(step)
      else :
         raise ValueError( "the \"limit\" list is not of appropriate length" )

      # axis is a list (of same size as _step) where each element is a np array
      # containing the nodes where "data" values are set
      for i in range(self._rank) :
         self.axis.append(np.linspace(limit[i][0], limit[i][1], self._size[i]))



   def _extractData(self,
                    data,
                    limit,
                    chunk,
                    stride) :

      """
      this method chunk & stride the data and correspondingly the axis
      both data & axis are modified

      input : limit is the list with the limit of the input data

              chunk is a list for min & max values (in physical units)
              between which data are to be kept

              stride is a list of integer containing the stride
              for each direction. if None, all data are kept

      """

      # manage iChunkMin & iChunkMax (containing index) from chunk
      if chunk == None :
         iChunkMin = []
         iChunkMax = []
         for i in range(self._rank) :
            iChunkMin.append(int(limit[i][0]/self._step[i]))
            iChunkMax.append(int(limit[i][1]/self._step[i]))

      else :
         if self._rank != chunk.__len__() :
            raise ValueError( "the \"chunk\" list is not of appropriate length" )

         else :
            iChunkMin = [None]*self._rank
            iChunkMax = [None]*self._rank
            for i in range(self._rank) :
               if chunk[i] == None :
                  iChunkMin[i] = int(limit[i][0]/self._step[i])
                  iChunkMax[i] = int(limit[i][1]/self._step[i])
               else :
                  if chunk[i][0] == None :
                     iChunkMin[i] = int(limit[i][0]/self._step[i])

                  else :
                     iChunkMin[i] = int(chunk[i][0]/self._step[i])

                  if chunk[i][1] == None :
                     iChunkMax[i] = int(limit[i][1]/self._step[i])

                  else :
                     iChunkMax[i] = int(chunk[i][1]/self._step[i])

      # manage stride
      if stride == None :
         stride = [1]*self._rank

      else :
         # replace all occurence of None by 1 in stride
         stride = [1 if x == None else x for x in stride]

      # make a list with slices along each axis specifying
      # min, max & stride
      obj = []
      for i in range(self._rank) :
         mySlice = slice(iChunkMin[i], iChunkMax[i]+1, stride[i])
         obj.append(mySlice)
         self.axis[i] = self.axis[i][mySlice]

      self.data = data[tuple(obj)]



   def _sliceData(self,
                  section) :

      """
      this method slice a given data in a/several given direction

      input : data is a numpy array of any rank

              section is a list of length "rank". if None, no slicing
              is performed ; if float value (in physical unit) a slice
              will be performed at the index associated to the float value

      output : a numpy array of rank smaller or equal as the one of data

      """


      if section == None :
         pass

      else :
         if self._rank != section.__len__() :
            raise ValueError( "the \"section\" list is not of appropriate length" )

         else :
            # build isection (index) from section (physical unit)
            isection = []
            for i, val in enumerate(section) :
               if val == None :
                  isection.append(None)
               else :
                  isection.append(int(val/self._step[i]))

            # list of axis to accordingly remove
            iaxis2remove = []

            # iterate along each axis, & slice if needed
            for i, val in enumerate(isection) :
               if val != None :
                  self.data = np.take(self.data, val, axis=i)
                  iaxis2remove.append(i)

            # this list has to be reversed because once an argument is poped
            # the index of the remaining values are modified... so values
            # need to be poped from last to first
            for i, val in reversed(list(enumerate(iaxis2remove))) :
               self.axis.pop(val)
               self._step.pop(val)


      self._rank = self.data.ndim
      self._size = self.data.shape
      #self._step = []

