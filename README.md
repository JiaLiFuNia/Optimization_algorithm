# 《最优化方法及其Matlab程序设计》与《微分方程数值方法》相关算法（Matlab）
> 
> 参考书目：
> 
>《最优化方法及其Matlab程序设计》马昌风 科学出版社 
> 
>《微分方程数值方法（第二版）》胡建伟、汤怀民 科学出版社 

## 《最优化方法及其Matlab程序设计》
### 1、 精确与非精确线搜索方法
- #### 进退法 [SF2_1](./YC/0/SF2_1.m)、[XT2_3](./YC/0/XT2_3.m)、[XT2_4](./YC/0/XT2_4.m)、[XT2_5](./YC/0/XT2_5.m)
- #### 抛物线法 [XT2_3](./YC/0/XT2_3.m)、[XT2_4](./YC/0/XT2_4.m)、[XT2_5](./YC/0/XT2_5.m)
- #### 0.618法 [XT2_1](./YC/0/XT2_1.m)、[XT2_2](./YC/0/XT2_2.m)、[XT2_5](./YC/0/XT2_5.m)
- #### Wolfe准则
- #### Armijo准则 [XT2_6](./YC/0/XT2_6.m)、[XT3_4](./YC/1/XT3_4.m)、[XT3_5](./YC/1/XT3_5.m)、[XT3_6](./YC/1/XT3_6.m)、[XT3_6_1](./YC/1/XT3_6_1.m)、[XT3_7](./YC/1/XT3_7.m)

### 2 、牛顿法与最速下降法
- #### 牛顿法 [SF3_2](./YC/1/SF3_2.m)、[其余](./YC/1/)
- #### 最速下降法 [SF3_1](./YC/1/SF3_1.m)、[其余](./YC/1/)
- #### 牛顿-最速下降法混合算法
- #### 修正牛顿法
- #### FR共轭梯度法
- #### 拟牛顿法
  
***

## 《微分方程数值方法》
### 1、 常微分方程初值问题的数值解法
- #### Euler法 [Euler_Method](./WF/Euler_Trapezoidal_Method.m)
- #### 梯形法 [Trapezoidal_Method](./WF/Euler_Trapezoidal_Method.m)
- #### 三阶龙格库塔法 [Runge Kutter](./WF/RungeKutter_Adams_PECE_PMCEMC.m)
- #### 四阶龙格库塔法 [Runge Kutter](./WF/RungeKutter_Adams_PECE_PMCEMC.m)
- #### Adams（显示格式） [Adams](./WF/RungeKutter_Adams_PECE_PMCEMC.m)
- #### PECE 模式 [PECE](./WF/RungeKutter_Adams_PECE_PMCEMC.m)
- #### PMECME 模式 [PMECME](./WF/RungeKutter_Adams_PECE_PMCEMC.m)