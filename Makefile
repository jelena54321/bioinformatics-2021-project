venv: venv/bin/activate

venv/bin/activate: requirements.txt
	python3 -m venv venv
	. ./venv/bin/activate && pip install pip --upgrade
	. ./venv/bin/activate && pip install -r requirements.txt

install: venv

clean:
	rm -rf venv
