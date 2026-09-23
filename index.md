## CCC for repeated measures

**rm_ccc_macro_v4.sas** is a SAS macro for estimating the concordance correlation coefficient (CCC) for repeated and non-repeated measurements.

The macro implements two approaches for estimating the CCC:

- a variance-components approach based on a linear mixed-effects model;
- a U-statistic approach for appropriate balanced designs.

### Version 4

Version 4 includes corrections and improvements to the variance-components implementation and allows the user to specify whether the observer-by-subject interaction is included in the model.

The `interaction` parameter controls the inclusion of the subject-by-observer interaction:

- `interaction=1`: includes the interaction term (default).
- `interaction=0`: excludes the interaction term, which may be useful when the interaction variance is negligible or its estimation causes numerical instability.

### Download

[Download rm_ccc_macro_v4.sas](https://github.com/jlcarrascoUB/rm_ccc_macro/raw/refs/heads/main/rm_ccc_macro_v4.sas)

[GitHub repository](https://github.com/jlcarrascoUB/rm_ccc_macro)

### Reference

Carrasco JL, Phillips BR, Puig-Martinez J, King TS, Chinchilli VM.  
**Estimation of the concordance correlation coefficient for repeated measures using SAS and R.**  
*Computer Methods and Programs in Biomedicine*. 2013;109(3):293–304.
