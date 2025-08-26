# protos/cli/cli.py
import typer
from .init_data import init_app  # Import sub-app from init_data.py
from .clean import clean_command  # Import individual command function from clean.py
from .register import register_app  # Import sub-app from register.py

app = typer.Typer(name="protos", help="CLI for the protos library.")

# Add an individual command directly (e.g., version)
@app.command()
def version():
    """Show the version."""
    typer.echo("protos version 1.0")  # Replace with actual version, e.g., from __version__

# Add individual imported commands
app.command()(clean_command)  # Adds the clean command

# Add sub-apps for grouped commands
app.add_typer(init_app, name="init")  # Makes commands like 'protos init ...'
app.add_typer(register_app, name="register")  # Makes 'protos register ...'

if __name__ == "__main__":
    app()