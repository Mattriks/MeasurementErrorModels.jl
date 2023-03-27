

# Theory

### General EIV model

```math
\begin{equation}
\begin{aligned}
    &\vb*{z}_t~=~ \vb*{z}^*_t + \vb*{q}_t \\
    &\vb*{z}^*_t ~=~ \left[ y^*_t \quad \vb*{x}^*_t\right] \\
    &\vb*{q}_t ~=~ \left[v_t \quad  \vb*{a}_t \right] \\
    &y_t^* ~=~ b_0 + \vb*{x}^*_t\vb*{B}^* + \epsilon_t
\end{aligned}
\end{equation}
```


```math
\begin{align}
\begin{bmatrix} Z_{11} \\ Z_{21} \\ \vdots  \\ Z_{n1}\end{bmatrix} &=~ 
    \begin{bmatrix}
        1 & \vb* x_1 \\ 1 & \vb* x_2 \\ \vdots \\ 1 & \vb* x_n 
    \end{bmatrix}
    \begin{bmatrix} b_0 \\ \vb*{B}^*    \end{bmatrix} + 
    \begin{bmatrix}\mathbb{I}_n \quad -(\vb*{B}^*{'} \otimes \mathbb{I}_n)    \end{bmatrix}
    \begin{bmatrix} \vb*{V} \\ vec(\mathbf A)   \end{bmatrix}
\end{align}
```



### Maximum likelihood




```math
\begin{align}
\vb B_\text{wtls} &=~\left[ (\vb X-\vb A)'\vb W^{-1}\vb X \right] ^{-1}\left[(\vb X-\vb A)\vb W^{-1}\vb Y\right] \\   
\Sigma_{B}^\text{wtls} &=~\left[ (\vb X-\vb A)'\vb W^{-1}\vb X \right]^{-1}\mathbf{M}\left[ (\vb X-\vb A)'\vb W^{-1}\vb X \right]^{-1'}
\end{align}
```
where ``\mathbf{M}=(\vb X-\vb A)'\vb W^{-1}\Sigma \vb W^{-1}(\vb X-\vb A)``.



### Bayesian solution

```math
\begin{align}
    \vb B &~=~ \left({\Sigma_{B}^\text{wtls}}^{-1} ~+~ \Sigma_0^{-1}\right)^{-1}\left({\Sigma_{B}^\text{wtls}}^{-1}\vb B_\text{wtls}~+~\Sigma_0^{-1}\vb B_0\right) \\
    \Sigma_{B} &~=~\left( {\Sigma_{B}^{\text{wtls}}}^{-1}~+~\Sigma_0^{-1}\right)^{-1}
\end{align}
```

### Reference

Fischer, M (2023) Measurement Error Proxy System Models: MEPSM v0.2. link to preprint