# some things the BCTM can do
 
##  [Leukemia Survival](leukemia)
**Proportional Hazards with Spatial Frailties**

![image](leukemia/leuk_ph.png)


Based on a dataset on acute myeloid leukemia survival ([[1]](#1)).

* analyzing impact of prognostic factors *age*, *sex*, *white blood cell count*  and *Townsend score*(indicating less affluent residential areas for higher values), 
* investigate spatial patterns for 24 administrative regions in North West England; . * reference is the minimum extreme value distribution resulting in a proportional hazards model

---

##  [Lung Cancer Survival](veteran)
**Semiparametric (Non-)Proportional Odds with Censoring**

![image](veteran/vet_densities.png)

Based on a Veteran’s Administration lung cancer trial dataset ([[2]](#2)).

* analyzing odds of survival dependent on Karnofsky Performance Score and different lung cancer types
* partial right-censoring

---

## References
<a id="1">[1]</a>
Henderson, R., Shimakura, S. and Gorst, D. (2002).
Modeling spatial variation in leukemia survival data.
Journal of the American Statistical Association 97(460): 965-972.

<a id="2">[2]</a>
Prentice, R. L. (1973).
Exponential survivals with censoring and explanatory variables.
Biometrika 60(2): 279-288.