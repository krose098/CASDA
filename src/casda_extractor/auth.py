"""CASDA authentication helpers with keyring-backed credential storage."""

from __future__ import annotations

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


def _inject_into_astroquery_keyring(casda_obj: Any, username: str, password: str) -> bool:
    """
    Store *password* in whichever keyring slot astroquery's login() reads from.

    astroquery's BaseQuery.login() looks up passwords from the system keyring
    using ``self._server`` (or a similarly named attribute) as the service key.
    If we populate that slot before calling login(), astroquery finds the
    password there and skips its own interactive prompt.

    Returns True if the password was successfully stored.
    """
    if not _KEYRING_AVAILABLE:
        return False
    for attr in ("_login_url", "_server", "LOGIN_URL", "SERVER", "server"):
        svc = getattr(casda_obj, attr, None)
        if isinstance(svc, str) and svc:
            try:
                keyring.set_password(svc, username, password)
                return True
            except Exception:  # noqa: BLE001
                pass
    return False


def get_casda(username: str | None = None, password: str | None = None) -> Casda:
    """
    Return a logged-in Casda instance, reusing an existing session if possible.

    Credential resolution order:
    1. ``username`` / ``password`` arguments.
    2. System keyring (if available and populated).
    3. Interactive prompt — username only; password is handled by astroquery.

    When we already have both credentials (from args or our keyring), the
    password is injected into astroquery's own keyring slot so that
    ``casda.login()`` proceeds without prompting.  When we only have the
    username, astroquery prompts for the password itself (matching the
    behaviour of a plain notebook ``casda.login(username=username)`` call).

    Parameters
    ----------
    username : str, optional
        OPAL username (email). Falls back to keyring / prompt if not supplied.
    password : str, optional
        OPAL password. Falls back to keyring / silent injection if not supplied.

    Returns
    -------
    Casda
        Authenticated astroquery Casda object.
    """
    from astroquery.casda import Casda  # lazy import — not needed at module load time

    global _casda_instance  # noqa: PLW0603
    if _casda_instance is not None:
        return _casda_instance

    # 1. Try our keyring for stored credentials
    if username is None or password is None:
        kr_user, kr_pass = load_credentials()
        if kr_user and kr_pass:
            username = username or kr_user
            password = password or kr_pass
            print(f"Using stored credentials for '{username}'.")

    # 2. Prompt for username if still missing
    if username is None:
        username = input("Enter your OPAL username (email): ").strip()

    # 3. Build the Casda instance
    casda = Casda()

    # 4. If we already have the password, inject it into astroquery's keyring
    #    slot so login() can find it without asking the user again.
    #    If injection fails (keyring unavailable or wrong attribute name),
    #    astroquery will fall through to its own interactive password prompt —
    #    which is exactly what a plain notebook login() call does.
    if password is not None:
        _inject_into_astroquery_keyring(casda, username, password)

    # 5. Login — astroquery handles the password lookup / prompt internally.
    #    We intentionally do NOT pass `password` here: CasdaClass._login()
    #    does not accept that kwarg and raises TypeError if it receives one.
    casda.login(username=username)
    _casda_instance = casda

    return casda


def reset_casda() -> None:
    """Clear the cached Casda session (forces re-login on next call)."""
    global _casda_instance  # noqa: PLW0603
    _casda_instance = None
