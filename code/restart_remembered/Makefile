union.trees : union.py init.trees rerun.trees
	python union.py

reset.trees : reset.py init.trees
	python reset.py

init.trees : run.slim
	slim run.slim

rerun.trees : rerun.slim reset.trees
	slim rerun.slim
