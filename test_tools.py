import unittest
import tools
import nose.tools as nt
import numpy as np
import pandas as pd



class TestToolsFunctions():
    def test_solve_for_u(self):
        e = np.random.uniform(0, 1)
        T = np.random.uniform(3282.3503, 3.46896e13)
        tau = np.random.uniform(3282.3503, 3.46896e13)

        #when t=tau we expect u=0
        #numerical solver gives approximately 0
        temp = tools.solve_for_u(tau, tau, T, e)
        nt.assert_almost_equal(0, temp)

        e = 0
        t = np.random.uniform(0, T)
        temp = tools.solve_for_u(t, tau, T, e)
        #when e=0 the function is linear
        nt.assert_almost_equal(2*np.pi*(t-tau)/T, temp)
    
    def test_get_data(self):
        temp = tools.get_data('test_data.tbl', '')
        nt.assert_equal(temp[0], pd.Series([0, 86400]))


    def test_find_files_for_star(self):
        nt.assert_raises(ValueError, tools.find_files_for_star, 'CrazyFrog')

