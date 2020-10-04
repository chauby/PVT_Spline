
#-*-coding:UTF-8 -*-
#%% ============== Import packages
import numpy as np

class PVT_Spline():
    def __init__(self, name='None'):
        """
        :param: name
            the name of the PVT_spline object
        """
        super(self.__class__, self).__init__()
        self.name = name

    def PVT_splineTwoPoints(self, via_point_0, via_point_1):
        """
        :param: via_point_0: float
            the first via points
        :param: via_point_1: float
            the second via points
        """
        p0 = via_point_0[0]
        v0 = via_point_0[1]
        t0 = via_point_0[2]

        p1 = via_point_1[0]
        v1 = via_point_1[1]
        t1 = via_point_1[2]

        T = t1-t0
        d = p0
        c = v0
        b = (3*p1 - v1*T - 2*v0*T - 3*p0)/(T**2)
        a = (-2*p1 + v1*T + v0*T + 2*p0)/(T**3)
        PVT_param = [a, b, c, d]

        return PVT_param

    def PVT_spline(self, via_points):
        """
        :param: via_points: float, np.array
            a set of the np.array with size Nx3, where N is the length of via points
        """
        queue_len = via_points.shape[0]
        if queue_len < 2:
            assert('point num too small')

        PVT_param_queue = np.zeros((queue_len-1, 4))

        for i in range(0, queue_len-1):
            PVT_param_queue[i,:] = self.PVT_splineTwoPoints(via_points[i,:], via_points[i+1,:])
        
        return PVT_param_queue

    def PVT_getPVTValue(self, PVT_param_queue, via_points, t):
        """
        :param: PVT_param_queue: float, list
            a set of the list object
        :param: via_points: float, np.array
            a set of the n[].array with size Nx3, where N is the length of via points
        """
        t_start = via_points[0, 2]
        t_end = via_points[-1, 2]
        PVT_point = np.zeros((1, 3))

        if t <= t_start:
            PVT_point_tmp = [via_points[0,0], via_points[0, 1], t]
        elif t >= t_end:
            PVT_point_tmp = [via_points[-1,0], via_points[-1, 1], t]
        else:
            for i in range(1, via_points.shape[0]):
                if t < via_points[i, 2]:
                    a = PVT_param_queue[i-1, 0]
                    b = PVT_param_queue[i-1, 1]
                    c = PVT_param_queue[i-1, 2]
                    d = PVT_param_queue[i-1, 3]
                    t0 = via_points[i-1, 2]

                    p = a*(t-t0)**3 + b*(t-t0)**2 + c*(t-t0) + d
                    v = 3*a*(t-t0)**2 + 2*b*(t-t0) + c
                    PVT_point_tmp = [p, v, t]
                    break

        PVT_point[0,0] = PVT_point_tmp[0]
        PVT_point[0,1] = PVT_point_tmp[1]
        PVT_point[0,2] = PVT_point_tmp[2]
        return PVT_point


# %% test the PVT method
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    via_points = np.zeros((5,3)) # N x 3, N is the length of via_points，3 is the dimension

    via_points[0,0] = 0 # pos
    via_points[0,1] = 0 # velocity
    via_points[0,2] = 0 # time

    via_points[1,0] = 2
    via_points[1,1] = -0.8
    via_points[1,2] = 1

    via_points[2,0] = 3
    via_points[2,1] = 0
    via_points[2,2] = 2.5

    via_points[3,0] = 1
    via_points[3,1] = 1.2
    via_points[3,2] = 4.5

    via_points[4,0] = 2
    via_points[4,1] = 0
    via_points[4,2] = 5.5

    pvt_spline = PVT_Spline()
    pvt_param_queue = pvt_spline.PVT_spline(via_points)
    data_len = 100
    pvt_trajectory = np.zeros((data_len, 3))
    time = np.linspace(0, via_points[-1, 2], data_len)
    for i,t in enumerate(time):
        pvt_trajectory[i,:] = pvt_spline.PVT_getPVTValue(pvt_param_queue, via_points, t)

    plt.figure()
    plt.plot(pvt_trajectory[:,2], pvt_trajectory[:,0], label='PVT spline')
    plt.plot(via_points[:,2], via_points[:,0], 'ro', label='via point')
    plt.title('PVT spline demo')
    plt.legend()
    plt.grid()
    plt.show()

# %%
