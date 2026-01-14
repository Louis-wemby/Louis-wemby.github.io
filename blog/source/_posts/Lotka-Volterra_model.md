---
title: 动态系统稳定性与 Lotka-Volterra 模型浅略研究
date: 2026-01-13
updated: 2026-01-13
description: A brief discussion about Lotka-Volterra model.
top_img: https://i.natgeofe.com/n/ee1a256e-5a93-4ddd-891b-6e3fb341cd83/01predators-phenomena.jpg
cover: https://oaknationalacademy-res.cloudinary.com/image/upload/v1700515161/xvjpxj8bmayxkhjawqkk.png
categories:
  - Exploration
tags: 
  - Modeling
  - Dynamics
comments: true
---

## 稳定性理论（Stability theory）

### 简介

稳定性理论研究微分方程的解及动态系统的轨迹在初始条件有微小扰动时的稳定性.在某初始平衡条件（Equilibrium）下，给系统向某一方向施加一微小扰动，若系统最终可能偏离平衡状态而无法回到原始位置，则该平衡点是**不稳定**的.反之，若为稳定平衡点，则最终能回到原始位置，如摆线的振荡、阻尼振子等系统.

<figure>
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/19/DampedCosine.svg/500px-DampedCosine.svg.png" alt="阻尼振荡">
    <figcaption>衰减的阻尼正弦波.来自：Wikipedia.</figcaption>
</figure>

### 稳定性的判断

动态系统可用微分方程来描述.对于一维的常微分方程如 $\frac{dx}{dt} = f(x)$，其平衡点 ${x^*}$ 满足 $f(x^\ast) = 0$.若 $f'(x^\ast) < 0$，则平衡点是稳定的；若 $f'(x^\ast) > 0$，则平衡点是不稳定的.值得注意的是，若 $f'(x^\ast) = 0$（如下图第二个 stable 处），则需要其他方法判断稳定性.

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/mainhttps://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/1D_stability.png?raw=true" alt="一维稳定性判断">
    <figcaption>一维动态系统的稳定性判断.</figcaption>
</figure>

而对于二维系统 $\frac{dx}{dt} = f_1(x_1, x_2)$ 和 $\frac{dy}{dt} = f_2(x_1, x_2)$，其平衡点 $(x^\ast, y^\ast)$ 满足 $f_1(x^\ast, y^\ast) = 0$ 和 $f_2(x^\ast, y^\ast) = 0$.通过线性化方法，可以构造类似雅可比矩阵（Jacobian matrix）：

$$
\mathbf{J} - \lambda \mathbf{I} =
\begin{bmatrix}
\frac{\partial f_1}{\partial x_1} - \lambda &
\frac{\partial f_1}{\partial x_2} \\\\
\frac{\partial f_2}{\partial x_1} &
\frac{\partial f_2}{\partial x_2} - \lambda
\end{bmatrix}
\tag {1}
$$

通过求解特征值 $\lambda$，可以判断平衡点的稳定性：
- 若实部均为**负**，则平衡点是**稳定**的（Stable Node）；
- 若实部均为**正**，则平衡点是**不稳定**的（Unstable Node）；
- 若实部**异号**，则平衡点是**鞍点**（Saddle Point），亦不稳定.

即 $(x_1,x_2)$ 为稳定点的充要条件为：

当

$$
\text{Det} \mathbf{J} = 0 \tag{2}
$$

时，

$$
\text{Re}\lambda _{1,2} < 0 \tag{3}
$$

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/mainhttps://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/2D_stability.png?raw=true" alt="二维稳定性判断">
    <figcaption>二维动态系统的稳定性判断.</figcaption>
</figure>

### 李雅普诺夫稳定（Lyapunov stability）

李雅普诺夫稳定常用来描述一个动力系统的稳定性.它是指系统在初始条件有微小扰动时，解的轨迹仍然保持在某一邻域内.具体来说，对于每个 $\epsilon > 0$，均存在 $\delta > 0$，使得当初始条件满足 $\lVert x(0) - x^\ast \rVert < \delta$ 时，对于所有 $t \geq {t_0}$，都有 $\lVert x(t) - x^\ast \rVert < \epsilon$，则称平衡点 $x^\ast$ 是李雅普诺夫稳定的.

在李雅普诺夫稳定的基础上，还可以定义渐近稳定（Asymptotic stability）和指数稳定（Exponential stability）等更强的稳定性概念：
- 渐近稳定：若平衡点李雅普诺夫稳定，且系统的解随着时间趋向于平衡点，即 $\lim_{t \to \infty} x(t) = x^\ast$ ，则称该平衡点是渐近稳定的.也就是说，渐近稳定不但能维持在平衡点附近，而且最后能收敛到平衡点.
- 指数稳定：若平衡点渐近稳定，且存在正数 $\alpha$ 和 $\beta$，使得 $\lVert x(t) - x^\ast \rVert \leq {\alpha} e^{-\beta t} \lVert x(0) - x^\ast \rVert$，则称该平衡点是指数稳定的.指数稳定表明状态函数的收敛速度不慢于某个指数函数的递减速度.

## Lotka-Volterra 模型

### 模型概述

Lotka-Volterra 模型由 Alfred Lotka（1925年）和 Vito Volterra（1926年）提出.该模型运用数学框架（常微分方程）描述自然界生态系统中捕食者与猎物种群数量随时间的动态变化，揭示两者间相互作用关系，对预测种群波动和理解生态系统的稳定性具有重要意义.

<figure>
    <img src="https://media.geeksforgeeks.org/wp-content/uploads/20240430104750/Lotka-Volterra-Model-(1).png" alt="Lotka-Volterra 模型曲线示意图">
    <figcaption>Lotka-Volterra 模型曲线示意图.</figcaption>
</figure>

Lotka-Volterra 模型通常基于以下基本假设：
1. 没有捕食者的情况下，猎物数量呈线性或指数增长.
2. 猎物可捕获性是限制捕食者种群增长的唯一因素.
3. 没有年龄结构差异，能持续繁殖，捕食者与猎物间具有原子同质性.
4. 捕食者与猎物的相遇几率和捕食率直接相关.5.捕食者拥有与种群密度无关的恒定死亡率.

该模型可由一对一阶非线性常微分方程描述：

$$
\begin{cases}
\frac{dx}{dt} = bx - pxy \\\\ 
\frac{dy}{dt} = rxy - dy
\end{cases}
\tag{4}
$$

其中变量 $x$ 表示猎物的种群数量，变量 $y$ 表示捕食者的种群数量.参数 $b$ 表示没有捕食者的情况下猎物种群的自发增长率.参数 $d$ 表示没有猎物的情况下捕食者的自然死亡率.乘积 $pxy$ 与 $rxy$ 分别表示由于捕食者存在猎物数量的净减少率，以及由于猎物存在捕食者数量的净增长率，其数值大小与二者种群数量相关.

### 稳定性分析
结合上述方程组，可提出以下问题：猎物与捕食者的种群数量将随时间如何变化？是否存在平衡点与稳定点？如果不稳定，能否形成振荡？持续的振荡需要满足何种参数条件，与初始种群数量是否有关系？

#### 理论计算

不难发现，方程组有两个平衡点 $(0,0)$ 与 $(\frac{d}{r}, \frac{b}{p})$.前者 $(0,0)$ 为不具备实际意义的平凡解，即猎物与捕食者的种群数量均为零（不稳定的鞍点），且当某一方的数量出现微小扰动时，便立即呈指数增长，从而偏离该平衡点.考虑第二个解 $(\frac{d}{r}, \frac{b}{p})$，令

$$
\begin{bmatrix}
b-py & -px \\\\
ry & rx - d
\end{bmatrix} = 0
\tag{5}
$$

可解得

$$
\lambda _{1,2} = \pm i\sqrt{bd} \tag{6}
$$

两个解的实部均等于零，可能为稳定点，也可能不是，需进一步讨论.尽管一般的非线性常微分方程组无法求出近似的通解，但对于本模型而言，在平衡态附近的微小扰动可求出近似通解

$$
\begin{cases}
x(t) = \frac{d}{r} + \sqrt{c}\ p A_0 \cos(\sqrt{b d}\ t + \phi) \\\\
y(t) = \frac{b}{p} + \sqrt{b}\ r A_0 \sin(\sqrt{b d}\ t + \phi)
\end{cases}
\tag {7}
$$

其中 ${A_0}$ 是与 $x$ 和 $y$ 初值有关的常数.解的周期为

$$
T = \frac{2\pi}{\sqrt{bd}} \tag{8}
$$

由此可见，对于微小扰动，系统会在平衡点附近做周期性振荡，即 $(\frac{d}{r}, \frac{b}{p})$ 为振荡中心，且振荡的频率与参数 $b$ 和 $d$ 有关.然而，这种振荡并非渐近稳定，因为解不会收敛到平衡点，而是持续围绕平衡点振荡.故此捕食者-猎物动态系统不存在稳定点.

#### 数值模拟
用 MATLAB 软件对相关数值进行模拟，以验证结果正确性，并研究不同初始条件（$x$、$y$ 的初值）下种群数量的变化趋势.相关代码如下：

```matlab
% Example of Lotka-Volterra Model Simulation
% Parameters
b=1;
p=0.01;
r=0.01;
d=2;

% ODE Solver
PredatorPrey=@(t,z)[b*z(1)-p*z(1)*z(2);r*z(1)*z(2)-d*z(2)];
[t,z]=ode45(PredatorPrey,[0,50],[250,50]);

% Plotting Results
subplot(1,2,1);
plot(t,z(:,1),"r-",t,z(:,2),"b-");
legend('Prey','Predator')
xlabel('T');
ylabel('Number of Predator/Prey');
title('Curve of Lotka-Volterra Model');
grid on;
subplot(1,2,2);
plot(z(:,1),z(:,2))
title('Phase Space')
grid on;
```

##### 验证 $(\frac{d}{r}, \frac{b}{p})$ 的平衡点特性

设置参数 $b=1$，$d=2$，$p=r=0.01$，则平衡点为 $(200,100)$.

设初值 $$(x_0,y_0) = (200,100)$$，在无扰动的情况下，种群数量保持不变，符合平衡点的特性.（下列所有结果图中，左边为种群数量变化曲线，红色代表猎物，蓝色代表捕食者；右边为相平面轨迹图.）

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/result1.png?raw=true" alt="Lotka-Volterra 模型模拟结果1">
    <figcaption>初值设为平衡点 $(200,100)$ 时的模拟结果.</figcaption>
</figure>

任取初值如 $$(x_0,y_0) = (199,99)$$，可见种群数量围绕平衡点做周期性振荡.

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/result2.png?raw=true" alt="Lotka-Volterra 模型模拟结果2">
    <figcaption>初值设为 $(199,99)$ 时的模拟结果.</figcaption>
</figure>

再如初值 $$(x_0,y_0) = (10,20)$$，同样可见种群数量围绕平衡点做周期性振荡.

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/result3.png?raw=true" alt="Lotka-Volterra 模型模拟结果3">
    <figcaption>初值设为 $(10,20)$ 时的模拟结果.</figcaption>
</figure>

由此可见，在不同初始条件下，种群数量均围绕平衡点 $(200,100)$ 做周期性振荡，验证了理论计算结果的正确性.然而，当初始条件偏离平衡点较近时，振荡幅度较小，相平面轨迹近似椭圆形闭合轨道；当初始条件偏离平衡点较远时，振荡幅度较大，相平面轨迹与椭圆形或圆形相差甚远，这是因为在鞍点 $(0,0)$ 附近取初值时，其相空间行为表现出一定的鞍点特征.上述结果表明振荡的幅度与初始条件有关，但振荡的频率（周期）仅与参数 $b$ 和 $d$ 有关.

##### 固定参数条件下不同初值对曲线变化趋势的影响

在平衡点附近各自取初值 $$(x_0,y_0) = (202,98)$$，

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/result4.png?raw=true" alt="Lotka-Volterra 模型模拟结果4">
    <figcaption>初值设为 $(202,98)$ 时的模拟结果.</figcaption>
</figure>

此时，发现两曲线并没有交点，两者种群数量分别在各自平衡点附近做振幅较小的周期性振荡.

再分别取$$(x_0,y_0) = (250,50)$$ 和 $$(x_0,y_0) = (1000,500)$$，可见种群数量虽然围绕平衡点做周期性振荡，但振荡幅度差异显著.尤其当偏离平衡位置较远时，振荡幅度可变得极其剧烈，不稳定的性质愈加明显.

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/result5.png?raw=true" alt="Lotka-Volterra 模型模拟结果5">
    <figcaption>初值设为 $(250,50)$ 时的模拟结果.</figcaption>
</figure>

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/result6.png?raw=true" alt="Lotka-Volterra 模型模拟结果6">
    <figcaption>初值设为 $(1000,500)$ 时的模拟结果.</figcaption>
</figure>

改变常微分方程组的参数值，使其均变为原来的两倍，即 $b=2$，$d=4$，$p=r=0.02$，此时平衡点仍为 $(200,100)$.

设初值 $$(x_0,y_0) = (250,50)$$.通过与前面结果的对比，可以发现种群数量依然围绕平衡点做周期性振荡，波动幅度（波峰和波谷）不变，但振荡频率更快，周期变小为原来的一半.

<figure>
    <img src="https://github.com/Louis-wemby/Louis-wemby.github.io/blob/main/blog/source/figure/posts/Lotka-Volterra_model/result7.png?raw=true" alt="Lotka-Volterra 模型模拟结果7">
    <figcaption>参数变为原来两倍，初值设为 $(250,50)$ 时的模拟结果.</figcaption>
</figure>

### 结论

通过理论分析与数值模拟，可以得出以下结论：Lotka-Volterra 模型的平衡点并非稳定点，而是振荡中心，种群数量会围绕该点持续振荡.在参数固定即平衡点不变的情况下，种群的初始数量会影响其变化趋势.其规律为：初始数量越接近平衡点，波动幅度越小；越偏离平衡点，波动幅度越大.因此，种群数量初始值与平衡点的差异在一定程度上可反映该系统的稳健性，即在微小干扰下种群数量是否会出现较大波动.然而，振荡的周期性不受种群初始数量的影响，所有变化曲线均表现出固定的周期性，其大小仅与参数 $b$ 和 $d$ 有关.

## 总结

Lotka-Volterra 模型是研究生态系统中捕食者与猎物种间关系这类动态系统的有力工具，能帮助我们了解不同物种在自然界中如何共存.但 Lotka-Volterra 模型也有其局限性，因为自然界是极其复杂的系统，影响种群数量与种间关系的因素是多方面的，难以被单一模型全部囊括.因此，在该模型的基础上，需要不断进行修正，如引入新的变量参数、考虑更高阶的模型，以求尽可能地捕获真实世界的复杂性.

## 参考文献

1. [Lyapunov stability](https://en.wikipedia.org/wiki/Lyapunov_stability)
2. GeeksforGeeks. (2024, April 30). Lotka-Volterra Model of Predator-Prey Relationship. GeeksforGeeks. <https://www.geeksforgeeks.org/lotka-volterra-model-of-predator-prey-relationship/>
3. 闵靖涛. (2018). Lotka-Volterra 模型的性质及生态学意义. 北京大学地球与空间科学学院地球物理系.