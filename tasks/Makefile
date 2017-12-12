include ../Makefile.inc

help: default
	@echo 
	@echo "You can also build and push specific images by name with:"
	@echo "make <task_name>"
	@echo "make push_<task_name>"
	@echo
	@echo "List available packages with 'make list'"


TASK=template

template:
	@if [[ "$(TASK)" =~ ^[a-z_][a-z0-9_]*$$ ]]; then \
		echo "Creating task folder for:  $(TASK)" ; \
		${TOOLS}/template_generator.py -t $(TASK); \
	else echo "Invalid task name, folders must be in lowercase separated by underscores, beginning with a letter"; \
	fi;


.PHONY: help   template
