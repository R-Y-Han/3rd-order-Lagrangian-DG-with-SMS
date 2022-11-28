---
html:
    toc: true
---

# A high-order Lagrangian DG hydrodynamic method for quadratic cells using a subcell mesh stablization scheme

Xiaodong Liu, Nathaniel R. Morgan, Donald E. Burton

## �㷨

### ���Ʒ���

$$
\begin{align}
\frac{dm}{dt} & = 0\\
\rho \frac{d\nu}{dt} & = \nabla \cdot \textbf{u}\\
\rho \frac{d\textbf{u}}{dt} & = \nabla \cdot \sigma\\
\rho \frac{d\tau}{dt} & = \nabla \cdot (\sigma \textbf{u})
\end{align}
$$
����$\sigma = -p \textbf{I}$��Ӧ��������$\tau = e + k$���������ܣ�$e$,$k$�ֱ����������ܺ��������ܡ���������ʽд�ɣ�
$$
\begin{equation}
\rho \frac{d\mathbb{U}}{dt} = \nabla \cdot \mathbb{H}
\end{equation}
$$

�ʵ��˶�����Ϊ
$$
\begin{equation}
\frac{dx}{dt} = \textbf{u}
\end{equation}
$$

### ����ӳ��

�ı�������ʹ�þŵ��ֵ��
$$
\begin{equation}
\textbf{x} (\xi,t) = \sum_{p} x_p(t) b_p(\xi)
\end{equation}
$$
![3rdordermesh](pics/3rdordermesh.png)

### ��ɢ

ѡȡTaylor������$\mathbf{\Psi}$������$(\xi_{cm},\eta_{cm})$���ǵ�Ԫ�����ģ��������$\mathbb{U}$�ɱ�ʾΪ��
$$
\begin{equation}
\mathbb{U}(\xi,t) = \mathbf{\Psi}(\xi) \cdot \mathbb{U}^k(t)
\end{equation}
$$
����
$$
\begin{equation}
\mathbb{U}^k = (\mathbb{U}_{cm}, \frac{\partial \mathbb{U}}{\partial \xi}\big|_{cm},\frac{\partial \mathbb{U}}{\partial \eta}\big|_{cm}, \frac{\partial^2 \mathbb{U}}{\partial \xi^2}\big|_{cm},\frac{\partial^2 \mathbb{U}}{\partial \xi \partial \eta}\big|_{cm}, \frac{\partial^2 \mathbb{U}}{\partial \eta^2}\big|_{cm})
\end{equation}
$$

$\mathbb{U}_{cm} = \overline{\mathbb{U}}$��$\mathbb{U}$�ڵ�Ԫ�ϵ�����ƽ��������Taylor���������Զ���ɣ�
$$
\begin{equation}
\begin{align*}
\Psi_1 & = 1\\
\Psi_2 & = \xi - \xi_{cm}\\
\Psi_3 & = \eta - \eta_{cm}\\
\Psi_4 & = \frac{(\xi-\xi_{cm})^2}{2} - \frac{1}{m} \int_{\omega(t)} \rho \frac{(\xi-\xi_{cm})^2}{2} d\omega\\
\Psi_5 & = \frac{(\eta-\eta_{cm})^2}{2}  - \frac{1}{m} \int_{\omega(t)} \rho \frac{(\eta-\eta_{cm})^2}{2} d\omega\\
\Psi_6 & = (\xi-\xi_{cm})(\eta - \eta_{cm}) - \frac{1}{m} \int_{\omega(t)} \rho (\xi-\xi_{cm})(\eta - \eta_{cm}) d\omega\\
\end{align*}
\end{equation}
$$

### DG

�ڷ���$\rho \frac{d \mathbb{U}}{dt} = \nabla \cdot \mathbb{H}$����ͬʱ���Ի�����$\Psi_m$Ȼ����֣����ø��ֹ�ʽ�õ���
$$
\begin{equation}
\begin{align*}
\sum_n M_{mn} \frac{d\mathbb{U}_n}{dt} & = \int_{\omega(t)} \nabla \cdot (\Psi_m \mathbb{H}) d\omega - \int_{\omega(t)} \mathbb{H} \cdot (\nabla \mathbb{H}) d\omega\\
& = \oint_{\omega(t)} \Psi_m \mathbb{H} \cdot \textbf{n} da - \int_{\Omega} \mathbb{H} \cdot (J^{-T}\nabla_{\xi} \Psi_m) |J| d\Omega
\end{align*}
\end{equation}
$$
����$M = (M_{mn})$���������󣬺�ʱ���޹ء�

����ϵĻ����ڲο���Ԫ���У��߽��ϵĻ���������ռ���У����Ǳ߽������ͨ���ڵ㴦��ֵʵ�ֵģ�����ʱֻ��Ҫ�ڵ㴦�ĵ�ֵ��

ʱ�����ʹ��RK3�����跽��Ϊ$M \nu = \mathbb{R}$����RK3Ϊ��
$$
\begin{equation}
\begin{align*}
& \nu^{s1} = \nu^n + \Delta t (M^{-1})^n \mathbb{R}^n\\
& \nu^{s2} = \frac{3}{4} \nu^n + \frac{1}{4} \nu^{s1} + \frac{\Delta t}{4}(M^{-1})^{s1} \mathbb{R}^{s1}\\
& \nu^{n+1} = \frac{1}{3} \nu^n + \frac{2}{3} \nu^{s2} + \frac{2 \Delta t}{3} (M^{-1})^{s2} \mathbb{R}^{s2}
\end{align*}
\end{equation}
$$

## RS

����corner force
$$
\begin{equation}
\textbf{F}_i^* = a_i \sigma^*_c \textbf{n}_i = a_i \sigma_c \textbf{n}_i + \mu_i a_i (\textbf{u}_V^* - \textbf{u}_c)
\end{equation}
$$
����$a_i$����ĳ�ֱ߳���$\textbf{n}_i$�ǣ�����ռ��У���λ�ⷨ�������±�$c$����corner�����ڽڵ㴦���������ع���ֵ��$\mu_i$���迹��$\textbf{u}_V^*$�ǽڵ��ٶȡ�

$$
\begin{equation}
\mu_i = \rho_z (c_z + \beta \frac{\gamma + 1}{2} (\textbf{u}_V^* - \textbf{u}_c) \cdot \textbf{n}_i)
\end{equation}
$$
�±�$z$����zone���ǵ�Ԫ������ֵ��

![cornerforce](pics/cornerforce.png)

### SMS

����������������񻮷�SMS����ԭ���ľŵ��ֵ���ı������񻮷�Ϊ�ĸ����������������е�ѹ��$\tilde{p}_c$ȥ����ԭ��RS�е�ѹ���������񻮷����£�
![meshdiscretization](pics/meshdiscretization.png)
![meshdiscretization2](pics/meshdiscretization2.png)

*���Ķ������񻮷ֵ���ʶ�������ģ���ÿ��������ȥ����һ��ƽ���ܶ� $\overline{\rho}_s$��������ܶȺ��ô���������ֵ�����������ƽ���ܶ� $\overline{\rho}_{\nu s}$�Ĳ�����������������Σ�Ȼ��ͨ����������޸�RS���õ���ѹ��ֵ��*

$\overline{\rho}_s$�ļ���Ϊ��
$$
\begin{equation}
\overline{\rho}_s = \frac{m_s}{|\omega_s(t)|} = \frac{\int_{\Omega} (\rho^0 |J|^0)_s d\Omega}{\int_{\Omega}(|J|)_s d\Omega}
\end{equation}
$$
��������������Ǵ�ǿ�����غ�õ��ģ�Ȼ����������ֵ�ع����ܶȼ���������ƽ���ܶȣ�
$$
\begin{equation}
\overline{\rho}_{\nu s} = \frac{\int_{\Omega} (\frac{1}{\nu} |J|)_s d\Omega}{|\omega_s(t)|}
\end{equation}
$$

Ȼ�󽫽ڵ���ܶ�$\rho_c$��Ϊ��
$$
\begin{equation}
\begin{align*}
\tilde{\rho}_c & = \rho_c + \chi \delta \rho_s\\
\delta \rho_s & = \overline{\rho}_s - \overline{\rho}_{\nu s}
\end{align*}
\end{equation}
$$
����$\chi$��һ������������û������˵����ȡΪ1������RS�е�$p_c$�͸ĳ�$\tilde{p}_c = EOS(\tilde{\rho}_c,e)$��RS�����ಿ�ֲ��䡣


### DG��Gauss����

֮ǰ�����õ���
$$
\begin{equation}
\sum_n M_{mn} \frac{d\mathbb{U}_n}{dt}
= \oint_{\omega(t)} \Psi_m \mathbb{H} \cdot \textbf{n} da - \int_{\Omega} \mathbb{H} \cdot (J^{-T}\nabla_{\xi} \Psi_m) |J| d\Omega
\end{equation}
$$

��������֣�$M_{mn}$�����ֶ����һ���ֱ����Gauss�����ʽ��������Ҫ��ע�߽�Ļ��֡�

#### Simpson�����ʽ��߽����

�Ա���$\nu$�Ļ���Ϊ����
$$
\begin{equation}
\sum_n M_{mn} \frac{d\nu_n}{dt} = (\sum_{i\in \omega(t)} \Psi_{m_i} a_i \textbf{n}_i \cdot \textbf{u}_V^* ) - [\cdots]
\end{equation}
$$
Simpson��ʽ��ϵ����������$a_i$�**�����$i$��Ҫ���������е�������߽磺**
![integration](pics/integration.png)

�Ա�$V_3V_4$Ϊ�����������ϵĻ�������д�����ǣ�
$$
\begin{equation}
(a_2 \textbf{n}_2)_{V_3} \cdot \textbf{u}_{V_3}^* + \sum_{i=1}^2 (a_i \textbf{n}_i)_{V_7^r} \cdot \textbf{u}_{V_7}^* + \sum_{i=1}^2 (a_i \textbf{n}_i)_{V_7^l} \cdot \textbf{u}_{V_7}^* + (a_1 \textbf{n}_1)_{V_4} \cdot \textbf{u}_{V_4}^*
\end{equation}
$$

> ʵ���Ͼ�������ÿ���������Ͻ���ͨ�����֣�Ȼ���ټ���һ����Ϊ�������ͨ�����֡��м�һ����$V_9$��Ϊ�������ڲ�����Ϊ�����ڲ�û�м����������$V_9$��ͨ������Ϊ0��

$a_i \textbf{n}_i$�ļ��㣺
ͬ����$V_3V_4$��Ϊ��������
$$
\int_{V_3V_4} d \textbf{s} = \int_{V_3V_4} \textbf{n} ds = \frac{1}{6} s_3 \textbf{n}_3 + \frac{4}{6} s_7 \textbf{n}_7 + \frac{1}{6} s_4 \textbf{n}_4
$$
����$s_i$��ʾ�ڵ�$V_i$��$|J|$��ֵ�����Ա���ֱ��ȡ��
$$
\begin{align*}
(a_2 \textbf{n}_2)_{V_3} & = \frac{1}{6} s_3 \textbf{n}_3\\
(a_1 \textbf{n}_1)_{V_4} & = \frac{1}{6} s_4 \textbf{n}_4\\
(a_1 \textbf{n}_1)_{V_7^r} & = (a_2 \textbf{n}_2)_{V_7^l} = \frac{2}{6} s_7 \textbf{n}_7
\end{align*}
$$

������߽��ϣ�
$$
(a_1 \textbf{n}_1)_{V_7^l} = -(a_2 \textbf{n}_2)_{V_7^r} = \zeta s_{7t} \textbf{n}_{7t}
$$
����$s_{7t}$Ϊ��$V_9V_7$��$|J|$�ڵ�$V_7$��ֵ��$\zeta$ȡ��$\frac{4}{6}$

### ������

��

### ����
>
> + $t_0$ʱ��
>   + ȷ�����񻮷�
>   + ��������
>   + �����ֵ
>   + ������������
> + $t_n$ʱ��
>   + RK stage i, i= 1,2,3��
>     + ����ڵ��ٶ�$u_p^*$
>     + �������м�����
>     + ����DG�Ļ���ֵ������������

***

## ��ֵ����
