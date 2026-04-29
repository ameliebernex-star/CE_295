# CE_295

## Formulation

### Variables
* time step $k$, container index $i$
* States: capacity $Q_i$, stored energy $E_i$, temperature $T_i$
* Decisions: charging and discharging power setpoints $c_i, d_i$, cooling effort $u_i$

State vector:

$$\vec{x} = [E_i(k), Q_i(k), T_i(k)]^\top \quad \forall \quad i \in [1,2], k\in[0,N]$$

Decision vector:

$$\vec{u} = [c_i(k), d_i(k), u_i(k)]^\top \quad \forall \quad i \in [1,2], k \in [0, N-1]$$

### Parameters
* time interval $\Delta t$
* degradation rate $\alpha$ in MWh capacity lost per MWh throughput
* temperature increase $\beta$ per MWh throughput
* temperature decrease $\gamma$ per cooling effort · hour
* inverter efficiency $\eta$
* temperature degradation factor $\kappa$
* capacity opportunity cost $\rho$
* replacement marginal cost $\sigma$
* baseline operating power $a$
* cooling operating power $b$

### Exogenous inputs
* expected value of real-time market price of electricity $\mathbb{E}[\lambda(k)]$

### State dynamics
* $E_i(k+1)=E_i(k)+\eta c_i(k) \Delta t-\frac{1}{\eta}d_i(k)\Delta t$
* $Q_i(k+1) = Q_i(k) - \alpha\phi\big(c_i(k)+d_i(k)\big)\Delta t$
	* where $\phi = 1 + \kappa \max \big( 0, T(k) - T_\text{ref} \big)$
* $T_i(k+1) = T_i(k) + \beta\big(c_i(k)+d_i(k)\big)\Delta t-\gamma u_i(k)\Delta t$

### Objective function
$$\min_{\vec{x},\vec{u}} J = \sum_{i=1}^2 \Bigg[C_N+ \sum_{k=0}^{N-1}  C_\text{arb}(k) + C_\text{oper}(k) + C_\text{repl}(k) + C_\text{opp}(k) \Bigg]$$

where

* arbitrage cost $C_\text{arb} = \mathbb{E}[\lambda(k)] \big(c_i(k)-d_i(k)\big)\Delta t$
* operating cost $C_\text{oper} = \mathbb{E}[\lambda(k)] \big(a+bu(k)\big)\Delta t$
* replacement marginal cost $C_\text{repl} = \sigma\Big(\alpha\phi\big(c_i(k)+d_i(k)\big)\Delta t\Big)$
* capacity opportunity cost $C_\text{opp} = \rho\Big(\alpha\phi\big(c_i(k)+d_i(k)\big)\Delta t\Big)$
* terminal cost $C_N = \quad 0 \text{ for } E_i(N) \ge \text{SoC}_N\cdot Q_i(N), \quad \infty \text{ for } E_i(N) \lt \text{SoC}_N\cdot Q_i(N)$

### Initial conditions
* Starting capacity: $Q_i(0) = Q_\text{new} \cdot \text{SoH}_0$
* Starting stored energy: $E_i(0) = Q_i(0) \cdot \text{SoC}_0$
* Starting temperature: $T_i(0)=T_\text{ref}$

### Boundary conditions
* Capacity must be non-negative: $Q_i(k) \ge 0$
* Energy operating limits: $\text{SoC}\_\text{min}\cdot Q_i(k)\le E_i(k)\le\text{SoC}\_\text{max}\cdot Q_i(k)$
* Temperature operating limits: $T_\text{min} \le T_i(k) \le T_\text{max}$
* Charging and discharging power setpoints must be non-negative can't exceed hard cap: $0 \le c_i(k),d_i(k) \le P_\text{max}$
* Cooling effort bounds: $0 \le u_i(k) \le 1$

### Terminal conditions
Captured by terminal cost.
