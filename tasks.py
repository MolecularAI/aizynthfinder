from invoke import task


@task
def build_docs(context):
    context.run("aizynthcli -h > ./docs/cli_help.txt")
    context.run("sphinx-apidoc -o ./docs ./aizynthfinder")
    context.run("sphinx-build -M html ./docs ./docs/build")


@task
def full_tests(context):
    cmd = "pytest --black --mccabe " \
          "--cov aizynthfinder --cov-branch --cov-report html:coverage --cov-report xml " \
          "tests/"
    context.run(cmd)

@task
def run_mypy(context):
    context.run("mypy --ignore-missing-imports --show-error-codes aizynthfinder")


@task
def run_linting(context):
    print("Running mypy...")
    context.run("mypy --ignore-missing-imports --show-error-codes aizynthfinder")
    print("Running pylint...")
    context.run("pylint aizynthfinder")