"""CASDA authentication helpers with keyring-backed credential storage."""

from __future__ import annotations

import getpass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from astroquery.casda import Casda

try:
    import keyring
    import keyring.errors

    _KEYRING_AVAILABLE = True
except ImportError:
    _KEYRING_AVAILABLE = False

_KEYRING_SERVICE = "casda-cli"
_KEYRING_USER_KEY = "casda_username"

_casda_instance: Any = None


# ---------------------------------------------------------------------------
# Credential management
# ---------------------------------------------------------------------------


def save_credentials(username: str, password: str) -> None:
    """
    Persist OPAL credentials in the system keyring.

    Parameters
    ----------
    username : str
        OPAL username (email address).
    password : str
        OPAL password.

    Raises
    ------
    RuntimeError
        If the ``keyring`` package is not installed.
    """
    if not _KEYRING_AVAILABLE:
        raise RuntimeError("Install the 'keyring' package to save credentials.")
    keyring.set_password(_KEYRING_SERVICE, _KEYRING_USER_KEY, username)
    keyring.set_password(_KEYRING_SERVICE, username, password)
    print(f"Credentials saved to system keyring for '{username}'.")


def load_credentials() -> tuple[str | None, str | None]:
    """
    Load OPAL credentials from the system keyring.

    Returns
    -------
    (username, password)
        Both ``None`` if nothing is stored or keyring is unavailable.
    """
    if not _KEYRING_AVAILABLE:
        return None, None
    try:
        username = keyring.get_password(_KEYRING_SERVICE, _KEYRING_USER_KEY)
        if username is None:
            return None, None
        password = keyring.get_password(_KEYRING_SERVICE, username)
        return username, password
    except keyring.errors.KeyringError:
        return None, None


def delete_credentials() -> None:
    """
    Remove stored OPAL credentials from the system keyring.

    Raises
    ------
    RuntimeError
        If the ``keyring`` package is not installed.
    """
    if not _KEYRING_AVAILABLE:
        raise RuntimeError("Install the 'keyring' package to manage stored credentials.")
    try:
        username = keyring.get_password(_KEYRING_SERVICE, _KEYRING_USER_KEY)
        if username:
            keyring.delete_password(_KEYRING_SERVICE, username)
        keyring.delete_password(_KEYRING_SERVICE, _KEYRING_USER_KEY)
        print("Stored credentials removed from system keyring.")
    except keyring.errors.KeyringError as exc:
        print(f"Could not remove credentials: {exc}")


# ---------------------------------------------------------------------------
# Session management
# ---------------------------------------------------------------------------


def get_casda(username: str | None = None, password: str | None = None) -> Casda:
    """
    Return a logged-in Casda instance, reusing an existing session if possible.

    Credential resolution order:
    1. ``username`` / ``password`` arguments.
    2. System keyring (if available and populated).
    3. Interactive prompt (username then password).

    If credentials were entered interactively and the keyring is available,
    the user is offered the option to save them for future runs.

    Parameters
    ----------
    username : str, optional
        OPAL username (email). Falls back to keyring / prompt if not supplied.
    password : str, optional
        OPAL password. Falls back to keyring / prompt if not supplied.

    Returns
    -------
    Casda
        Authenticated astroquery Casda object.
    """
    from astroquery.casda import Casda  # lazy import — not needed at module load time

    global _casda_instance  # noqa: PLW0603
    if _casda_instance is not None:
        return _casda_instance

    entered_manually = False

    # 1. Try keyring
    if username is None or password is None:
        kr_user, kr_pass = load_credentials()
        if kr_user and kr_pass:
            username = username or kr_user
            password = password or kr_pass
            print(f"Using stored credentials for '{username}'.")

    # 2. Interactive fallback
    if username is None:
        username = input("Enter your OPAL username (email): ").strip()
        entered_manually = True
    if password is None:
        password = getpass.getpass(f"Enter OPAL password for {username}: ")
        entered_manually = True

    casda = Casda()
    casda.login(username=username, password=password)
    _casda_instance = casda

    # Offer to save if entered manually and keyring is available
    if entered_manually and _KEYRING_AVAILABLE:
        try:
            save_now = input("Save credentials to system keyring for future runs? [y/N] ").strip().lower()
            if save_now == "y":
                save_credentials(username, password)
        except EOFError:
            pass  # Non-interactive environment — skip the prompt silently

    return casda


def reset_casda() -> None:
    """Clear the cached Casda session (forces re-login on next call)."""
    global _casda_instance  # noqa: PLW0603
    _casda_instance = None
