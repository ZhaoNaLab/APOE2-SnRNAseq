default.resources = list(h_vmem = '20G', queue = '7-day', email = 'YourEmail@gmail.com', conda.env = "../CondaEnv/EnvName", measure.memory = TRUE)
cluster.functions = makeClusterFunctionsSGE(template = "sge-simple.tmpl")

