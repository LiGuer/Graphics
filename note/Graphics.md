image is a matrix, and its element is vector ARGB

$$Img \in R^{H \times W \times C}, pix = \left(\begin{matrix}A \\ R \\ G \\ B\end{matrix}\right)$$

* Bsaic Geometry
  * Draw Line
    - Problem  
      (input) 2 points $\boldsymbol p_1, \boldsymbol p_2 \in \mathbb Z^2$  
      if $\boldsymbol p = \left(\begin{matrix}x \\ y\end{matrix}\right) \in \mathbb R^2$ in the line formed by $p_1$ and $p_2$, then

      $$(x_2 - x_1)(y - y_1) = (y_2 - y_1)(x - x_1)$$

      we need find all pixels $\boldsymbol p \in \mathbb Z^2$ that the line formed by $p_1$ and $p_2$ pass through, and minimize the error $\epsilon \ge 0$,

      $$(x_2 - x_1)(y - y_1) = (y_2 - y_1)(x - x_1) \pm \epsilon$$

      we set $\Delta_x, \Delta_y, a_x, a_y \in \mathbb Z$,
        $$\Delta = \left(\begin{matrix}\Delta_x \\ \Delta_y\end{matrix}\right) = \left(\begin{matrix}x_2 - x_1 \\ y_2 - y_1\end{matrix}\right)$$

      $$\begin{align*}
        \Rightarrow\quad 
        \min \quad& \epsilon = | a_y - a_x |  \\
        s.t. \quad
        & a_y = \Delta_x (y - y_1)  \\
        & a_x = \Delta_y (x - x_1)
      \end{align*}$$

      For higher dimensions, $\cdot^{(i)}$ refers to the value of the $i$-th dimension, the target is minimize the defferece $\epsilon$ among $a^{(i)}, i = 1,...,\dim$.

      $$\begin{align*}
      \min \quad& \epsilon = \sum_{i,j = 1:\dim, i \neq j} | a^{(i)} - a^{(j)} | \\
      s.t. \quad & a^{(i)} = \left(x^{(i)} - x_1^{(i)}\right) \cdot \prod_{d = \{1:\dim\}-\{i\}} \Delta^{(i)}
      \end{align*}$$


    - Algorithm: Bresenham Algorithm  
      The target of this algorithm is to solve the problem without using any real number.
      then, we find
      $$x \gets x \pm 1 \quad\Rightarrow\quad a_x \gets a_x \pm \Delta_y$$
      $$y \gets y \pm 1 \quad\Rightarrow\quad a_y \gets a_y \pm \Delta_x$$

      We propose the algorithm: ```while``` $x \gets x \pm 1$ and $a \gets a \pm \Delta_y$, ```if``` $a \% \Delta_x$ has increased by $1$ compared with the previous moment, ```then``` $y \gets y \pm 1$. (The selection of $\pm$ depend on the relative size of $x_1, x_2$ and $y_1, y_2$)
      
      Meanwhile, We set the dimension $d$ with the max $Δ_d$ as the the above $x$ to maintain continuity of the line and avoid impact of $Δ_d = 0$.

      ```cpp
      while x = x + 1
        a_ = a
        a = a + Δ_y

        if (a % Δ_x) - (a_ % Δ_x) >= 1 
          y = y + 1
      ```

      For higher dimensions, the algorithm is a similar as above.
      
  * Draw Circle, Sphere, Hyper-Sphere
    - Problem
      $$\begin{align*}
        x^2 + y^2 &= r^2  \tag{Circle, 2D}  \\
        x^2 + y^2 + z^2 &= r^2  \tag{Sphere, 3D} \\
        \|\boldsymbol p\|_2^2 = \sum_{i=1:\dim} p^{(i)2} &= r^2    \tag{any dimonsion}
      \end{align*}$$

      The target is to find all pixels $\boldsymbol p \in \mathbb Z^{\dim}$ that the boundary of Circle, Sphere or Hyper-Sphere pass through, and minimize the error $\epsilon \ge 0$. $x, y, r \in \mathbb Z$.
      $$\begin{align*}
        \min \quad& \epsilon = |x^2 + y^2 - r^2|  \\
        \min \quad& \epsilon = \left|\|\boldsymbol p\|_2^2 - r^2 \right|
      \end{align*}$$

    - Algorithm: Bresenham Algorithm   
      We only need to draw $\frac{1}{8}$ of a circle, and then we will get the whole circle through symmetry. If we start in $(0, r)$, then we draw in order  
      $$\left(\theta: \frac{\pi}{2} \to \frac{\pi}{4}, x: 0 \to \frac{r}{\sqrt{2}}, y: r \to \frac{r}{\sqrt{2}}, x \le y\right)$$  

      Therefore, the next position of each step is only $(x+1, y)$ and $(x+1, y-1)$. And, we only need to select the position with smaller error $\arg\min(\epsilon_1, \epsilon_2)$ by the positive and negative of $\delta$.

      $$\begin{align*}
        (x+1)^2 + y^2 - r^2 &= + \epsilon_1 \ge 0\\
        (x+1)^2 + (y-1)^2 - r^2 &= - \epsilon_2  \le 0
      \end{align*}$$

      $$\begin{align*}
        \delta &= + \epsilon_1 - \epsilon_2  \\
        &= ((x+1)^2 + y^2 - r^2) + ((x+1)^2 + (y-1)^2 - r^2)
      \end{align*}$$

      $$y \gets \left\{\begin{matrix}
          y \quad&; \delta_{x} < 0  \\
          y + 1 \quad&; \delta_{x} \ge 0
        \end{matrix}\right.$$

      Recursion formula and the initial value of $\delta$,

      $$\begin{align*}
        \delta_{(0,r)} &= 3 - 2 r  \\
        \delta_{x+1} &= \delta_{x} + \Delta \delta
      \end{align*}$$

      $$\begin{align*}
        \Delta \delta &= \delta_{x+1} - \delta_{x}  \\
        &= \left\{\begin{matrix}
          4 x + 6 \quad&; \delta_{x} < 0  \\
          4 x - 4 y + 10 \quad&; \delta_{x} \ge 0
        \end{matrix}\right.
      \end{align*}$$

      For higher dimensions, we only need to draw $\frac{1}{N}$ of a circle, through the positivity and negativity of axis $\frac{1}{2^{\dim}}$ and permutation of axis $\frac{1}{\dim !}$. And this area is divided by hyperplane $p^{(i)} = 0, p^{(i)} = p^{(j)}, \forall i, j \in \{1:\dim\}$.

      $$\frac{1}{N} = \frac{1}{2^{\dim} \cdot \dim !}$$

      We start from the initial point $(0, ..., 0, r)$, and satisfy:
      $$\begin{align*}
        \Delta p^{(i)} &\in \{0, +1\} \quad; \forall i \in \{1:\dim-1\}\\
        \Delta \boldsymbol p^{(1:\dim-1)} &\neq \boldsymbol 0  \\
        \Delta p^{(\dim)} &\in \{0, -1\}
      \end{align*}$$

      $$p^{(1)} \le p^{(2)} \le ... \le p^{(\dim)}$$

      Then, we search all feasible points satisfying $\epsilon \le r$ by ```queue``` and breadth-first search.

      $$\begin{align*}
      \Delta d &= \sqrt{\|p\|_2^2} - r\\
      &\in \left[-\frac{\sqrt{\dim}}{2}, \frac{\sqrt{\dim}}{2}\right] \\
      \Rightarrow\quad \epsilon &= \left|\|p\|_2^2 - r^2 \right|  \\
      &= \left|(\Delta d + r)^2 - r^2\right|  \\
      &= \left|\Delta d^2 + 2 r \Delta d + r^2 - r^2\right|\\
      &= \left|\Delta d^2 + 2 r \Delta d\right|  \\
      &\in [0,r \sqrt{\dim}  +\frac{\dim}{4}] \\
      \end{align*}$$

  * Draw Analytic Geometry of Arbitrary form

    We search all feasible points satisfying $f(\boldsymbol p) \le 0$ from a initial feasible point $\boldsymbol p_0$ by ```queue``` and breadth-first search.

  * Draw Triangle
    - Algorithm
      - Scanline 