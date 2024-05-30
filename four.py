import os
import scipy.constants as constants
import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt


#Чтение из файла
import urllib.request
import xml.etree.ElementTree as ET
url = 'https://jenyay.net/uploads/Student/Modelling/task_rcs.xml'
with urllib.request.urlopen(url) as response:
    dataread = response.read()
tree = ET.fromstring(dataread)
for elem in tree.iter('variant'):
    if elem.attrib['number'] == '2':
        D = float(elem.attrib['D'])
        fmin = float(elem.attrib['fmin'].replace('e9', '')) * 1e9
        fmax = float(elem.attrib['fmax'].replace('e9', '')) * 1e9

freq= range(int(fmin) , int(fmax), 10000000)

class RCS:
    def __init__(self, D, freq):
        self.r = D/2
        self.wave_length = 0
        self.k = 0
        self._freq_range_ = freq
    def a_n(self, n):
        return np.longdouble(special.spherical_jn(n, self.k * self.r))/self.h_n(n, self.k * self.r)

    def b_n(self, n):
        numerator = self.k * self.r * np.longdouble(special.spherical_jn(n - 1, self.k * self.r)) - n * np.longdouble(special.spherical_jn(n, self.k * self.r))
        divider = self.k * self.r * self.h_n(n - 1, self.k * self.r) - n * self.h_n(n, self.k * self.r)
        return np.divide(numerator, divider)
    
    def h_n(self, n, arg):
        return np.clongdouble(special.spherical_jn(n, arg) + 1j*special.spherical_yn(n, arg))
    
    def calculRCS(self):
        coef = self.wave_length**2 / constants.pi
        partForml = 0
        # оператор суммы в формуле c верхним пределом 50
        for n in range(1, 50):
            partForml += (-1)**n * (n + 0.5) * (self.b_n(n) - self.a_n(n))
        result = coef * abs(partForml) ** 2
        return result
    
    def GetData(self):
        self.data_freq = []
        self.data_lamda = []
        self.data_rcs = []
        for freq in self._freq_range_:
            # обновляем длину волны и волновое число
            self.wave_length = np.longdouble(constants.speed_of_light / freq)
            self.k = np.longdouble(2 * constants.pi / self.wave_length)
            # получаем значение ЭПР для новых параметров
            temp_rcs = self.calculRCS()
            self.data_freq.append(float(freq))
            self.data_lamda.append(float(self.wave_length))
            self.data_rcs.append(float(temp_rcs))
        return self.data_freq, self.data_lamda, self.data_rcs

class Output:
    def __init__(self, data_freq, data_lamda, data_rcs):
        self.data_freq= data_freq
        self.data_lamda= data_lamda
        self.data_rcs= data_rcs
        
    def save_to_txt(self, directory_name):
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
        file_path = os.path.join(directory_name, 'results.txt')
        with open (file_path, 'w') as file:
            for freq, lamda, rcs in zip (self.data_freq, self.data_lamda, self.data_rcs):
                line=f'{freq}\t{lamda}\t{rcs}\n'
                file.write(line)

    def plot_data(self):
        plt.plot(self.data_freq, self.data_rcs)
        plt.xlabel('Frequency')
        plt.ylabel('RCS')
        plt.title('RCS from frequency')
        plt.grid()
        plt.show()
        
calculator = RCS(D, freq)
freq,lamda,rcs = calculator.GetData()
output = Output(freq, lamda,rcs)
output.save_to_txt('results')
output.plot_data()

