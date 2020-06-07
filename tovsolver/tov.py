from scipy.interpolate import interp1d
from scipy.integrate import odeint
from scipy.optimize import root
import matplotlib.pyplot as plt
import pkg_resources

import numpy as np

from .constants import *

class TOV:
  def __init__(self, en_arr, p_arr, N=1000, add_crust=True, plot_eos=False):

    en_arr *= MeV_fm3_to_pa_cgs / c**2
    p_arr  *= MeV_fm3_to_pa_cgs

    sort_ind = np.argsort(p_arr)
    self.en_dens = interp1d(p_arr[sort_ind], en_arr[sort_ind], kind='cubic')

    sort_ind = np.argsort(en_arr)
    self.press = interp1d(en_arr[sort_ind], p_arr[sort_ind], kind='cubic')

    self.__en_arr = en_arr
    self.__p_arr = p_arr

    self.min_dens = np.min(en_arr)
    self.max_dens = np.max(en_arr)

    self.min_p = np.min(p_arr)
    self.max_p = np.max(p_arr)

    self.N  = N

    if add_crust:
      if(plot_eos):
        plt.plot(self.__en_arr / (MeV_fm3_to_pa_cgs / c**2),
                 self.__p_arr/MeV_fm3_to_pa_cgs ,
                 linestyle='-', label='original EOS')
      self.add_crust()
    if(plot_eos):
      plt.plot(self.__en_arr / (MeV_fm3_to_pa_cgs / c**2), 
               self.__p_arr/MeV_fm3_to_pa_cgs,
               linestyle='--', label='EOS with crust')

      plt.xscale('log')
      plt.yscale('log')
      plt.xlabel(r'${\rm \varepsilon~(MeV/fm^{3}) }$')
      plt.ylabel(r'${\rm P~(MeV/fm^{3}) }$')
      plt.legend()
      plt.show()

  def add_crust(self):
    crust_loc = pkg_resources.resource_filename(__name__, 'data/')
    # dir_name = os.path.dirname(__file__)
    
    baym_eos = np.genfromtxt(crust_loc + "Baym_eos.dat", 
                         dtype=float, skip_header=1, unpack=True, 
                         names=["en", "p", "nB",])

    P_crust = interp1d(baym_eos["en"], baym_eos["p"], kind = 'cubic')

    def eq_glue(n):
      return P_crust(n) - self.press(n)

    g = root(eq_glue, [44.*(MeV_fm3_to_pa_cgs / c**2)], options = {'maxfev' : 200})

    n_glue = g['x'][0]
    
    en_arr = []
    p_arr = []

    for i in range(len(baym_eos["p"])):
        if baym_eos["en"][i] < n_glue:
            en_arr.append(baym_eos["en"][i])
            p_arr.append(baym_eos["p"][i])
        else:
            break

    glue_ind = i

    for i in range(len(self.__p_arr)):
        if self.__en_arr[i] >= n_glue:
            en_arr.append(self.__en_arr[i])
            p_arr.append(self.__p_arr[i])

    en_arr = np.array(en_arr)
    p_arr = np.array(p_arr)

    self.min_dens = min(en_arr)
    self.min_p = min(p_arr)

    self.en_dens = interp1d(p_arr, en_arr, kind='cubic')
    self.press = interp1d(en_arr, p_arr, kind='cubic')

    self.__en_arr = en_arr
    self.__p_arr  = p_arr

    return

  def dedp(self, r, R_dep):
    e_R, p_R, m_R = R_dep

    p = p_R(r)
    dp = p * 0.005

    el_3 = self.en_dens(p - 3 * dp)
    el_2 = self.en_dens(p - 2 * dp)
    el_1 = self.en_dens(p - 1 * dp)
    er_3 = self.en_dens(p + 3 * dp)
    er_2 = self.en_dens(p + 2 * dp)
    er_1 = self.en_dens(p + 1 * dp)
    de_dp = (-1 / 60 * el_3 + 3 / 20 * el_2 - 3 / 4 * el_1 + 3 / 4 * er_1 - 3 / 20 * er_2 + 1 / 60 * er_3) / dp

    return de_dp

  def love_eq(self, param, r, R_dep):
    beta, H = param
    e_R, p_R, m_R = R_dep

    try:
      dummy = p_R(r)
    except ValueError:
      return [100000, 100000]

    de_dp = self.dedp(r, R_dep)

    dbetadr = H * (-2 * pi * G / c ** 2 * (
        5 * e_R(r) + 9 * p_R(r) / c ** 2 + de_dp * c ** 2 * (e_R(r) + p_R(r) / c ** 2)) \
                   + 3 / r ** 2 \
                   + 2 * (1 - 2 * m_R(r) / r * km_to_mSun) ** (-1) * (
                       m_R(r) / r ** 2 * km_to_mSun + G / c ** 4 * 4 * pi * r * p_R(r)) ** 2) \
              + beta / r * (
                  -1 + m_R(r) / r * km_to_mSun + 2 * pi * r ** 2 * G / c ** 2 * (e_R(r) - p_R(r) / c ** 2))
    dbetadr *= 2 * (1 - 2 * m_R(r) / r * km_to_mSun) ** (-1)

    dHdr = beta
    return [dbetadr, dHdr]

  def tov_eq(self, y, r):
    P, m = y

    if P < self.min_p or P > self.max_p:
      return [0., 0.]

    eden = self.en_dens(P)

    dPdr = -G * (eden + P / c ** 2) * (m + 4.0 * pi * r ** 3 * P / c ** 2)
    dPdr = dPdr / (r * (r - 2.0 * G * m / c ** 2))

    dmdr = 4.0 * pi * r ** 2 * eden

    return [dPdr, dmdr]


  def solve(self, c_dens, rmax=30e5, rtol=1.0e-5, dmrel=10e-12, dr=100):
    c_dens *= MeV_fm3_to_pa_cgs / c ** 2
    N = self.N
    r = np.arange(dr, rmax + dr, dr)

    if c_dens < self.min_dens or c_dens > self.max_dens:
      raise Exception('Central density: %8.4E is outside of the EoS range' %(c_dens))

    P = self.press(c_dens)
    eden = self.en_dens(P)
    m = 4.0 * pi * r[0] ** 3 * eden

    psol = odeint(self.tov_eq, [P, m], r, rtol=rtol)

    p_R, m_R = psol[:,0], psol[:,1]

    # find the boundary of the star by finding point
    # where the mass stops to increase

    diff = (m_R[1:] - m_R[:-1])/m_R[1:]
    ind = -1
    for i, dm in enumerate(diff):
      if dm < dmrel and m_R[i] != 0:
        ind = i
        break

    M = m_R[ind - 1]
    R = r[ind - 1]

    r   = r[:ind]
    p_R = p_R[:ind]
    m_R = m_R[:ind]
    
    e_R = self.en_dens(p_R)
    
    return R / 1e5, M / Msun, (r, e_R, p_R, m_R)

  def solve_tidal(self, c_dens, rmax=30e5, rtol=1.0e-4, dr=100):
    R, M, R_dep  = self.solve(c_dens, rmax, dr=dr)
    r, e_R, p_R, m_R = R_dep

    R *= 1e5
    M *= Msun

    e_R = interp1d(r, e_R, kind='cubic')
    p_R = interp1d(r, p_R, kind='cubic')
    m_R = interp1d(r, m_R, kind='cubic')

    beta0 = 2 * r[0]
    H0 = r[0] ** 2

    solution = odeint(self.love_eq, [beta0, H0], r, args=([e_R, p_R, m_R],), rtol=rtol)

    beta = solution[-1, 0]
    H = solution[-1, 1]

    y = R * beta / H

    C = compactness = M / R * km_to_mSun

    k2 = 8 / 5 * C ** 5 * (1 - 2 * C) ** 2 * (2 + 2 * C * (y - 1) - y) * (
          2 * C * (6 - 3 * y + 3 * C * (5 * y - 8)) + 4 * C ** 3 * (
            13 - 11 * y + C * (3 * y - 2) + 2 * C ** 2 * (1 + y)) + 3 * (1 - 2 * C) ** 2 * (2 - y + 2 * C * (y - 1)) * (
            np.log(1 - 2 * C))) ** (-1)

    return np.array([R / 1e5, M / Msun, C, k2, y, beta, H])




