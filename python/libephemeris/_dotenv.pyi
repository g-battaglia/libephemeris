from pathlib import Path
from typing import Any, Optional, Union, Dict, Literal
import os

def load_dotenv(path: Optional[Union[str, Path]]=None, *, override: bool=False) -> bool: ...
