from mcu.utils import elastic2D

# Ref: https://doi.org/10.1021/jacs.8b13075

# Elastic tensors
AlB6_1 = [379.9,438.4,23.1,159.6]
AlB6_2 = [383,375.3,33.6,132.1]
AlB6_3 = [395.1,401,44,173.6]
AlB6_4 = [229.9,194.2,7.1,80.3]
AlB6_5 = [242.2,171.1,15.6,57.1]
AlB6_6 = [149.3,92.2,24.9,63.2]

# define a list of elastic tensors
elastic_tensors = [AlB6_1,AlB6_2,AlB6_3,AlB6_4,AlB6_5,AlB6_6]

# Analyse
legend = ["AlB6_1", "AlB6_2", "AlB6_3", "AlB6_4", "AlB6_5", "AlB6_6",]
elastic2D.analyze(elastic_tensors)
elastic2D.plot_polar(elastic_tensors, young=True, legend=legend, figname='young', save=True)
elastic2D.plot_polar(elastic_tensors, young=False, legend=legend, figname='poisson', save=True)