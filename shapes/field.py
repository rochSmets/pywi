
import numpy as np
import shape


class Field(shape.Shape):

    """

    """

    def __init__(self,
                 limit,
                 data,
                 domain = None,
                 stride = None,
                 section=None,
                 shifts = None,
                 labels = None):

        limit = list(limit)

        super(Field, self).__init__(data, limit, domain, stride, section)

        # .. set the shifts for the axis
        self._setShifts(shifts)

        # .. set labels
        self._setLabels(section, labels)


    def _setShifts(self,
                   shifts) :


        if shifts == None :
            shifts = self.axis.__len__()*[0]

        else :
            for i, val in enumerate(shifts) :
                if val == None :
                    shifts[i] = 0.0
                else :
                    pass

        for i in range(self._rank) :
            self.axis[i] += shifts[i]


    def _setLabels(self,
                   section,
                   labels) :

        """
        should be virtual & defined in field, fourier... !!!!!!!!!!!!!!!!!!!!!!!!!
        """

        if labels != None :
            self.labels = labels

        else :
            self.labels = []
            s = ['x', 'y', 'z']

            if section == None :
                for i in range(self.data.ndim) :
                    self.labels.append('$'+s[i]+' / d_0$')

            else :
                for i in range(section.__len__()) :
                    if section[i] == None :
                        self.labels.append('$'+s[i]+' / d_0$')

                    else :
                        pass

