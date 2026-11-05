"""Standard logger integration and native binary size formatting."""

import inspect
import logging

import libephemeris as eph
from libephemeris import logging_config as config
import pytest


@pytest.fixture(autouse=True)
def restore_logging():
    logger = config.get_logger()
    previous = logger.level
    config.enable_logging(config.DEFAULT_LEVEL)
    yield
    config.set_log_level(previous)


def test_logger_is_the_standard_singleton_and_configuration_is_not_duplicated():
    first = config.get_logger()
    count = len(first.handlers)
    assert isinstance(first, logging.Logger)
    assert first is logging.getLogger("libephemeris")
    assert config.get_logger() is first
    assert len(first.handlers) == count
    assert first.name == config.LOGGER_NAME == "libephemeris"
    assert first.propagate is False
    assert config.LIBEPHEMERIS_LOG_LEVEL_ENV == "LIBEPHEMERIS_LOG_LEVEL"
    assert list(inspect.signature(config.set_log_level).parameters) == ["level"]
    for name in ("get_logger", "set_log_level", "disable_logging", "enable_logging", "format_file_size"):
        assert getattr(config, name) is getattr(eph, name)


def test_level_set_disable_and_default_enable():
    logger = config.get_logger()
    config.set_log_level(logging.DEBUG)
    assert logger.level == logging.DEBUG
    assert all(handler.level == logging.DEBUG for handler in logger.handlers)
    config.disable_logging()
    assert logger.level == logging.CRITICAL + 1
    assert not logger.isEnabledFor(logging.CRITICAL)
    config.enable_logging(logging.ERROR)
    assert logger.level == logging.ERROR
    config.enable_logging()
    assert logger.level == config.DEFAULT_LEVEL


def test_stdlib_logger_accepts_normal_handlers_and_message_arguments():
    records = []

    class Capture(logging.Handler):
        def emit(self, record):
            records.append(record)

    logger = config.get_logger()
    handler = Capture()
    logger.addHandler(handler)
    try:
        config.set_log_level(logging.DEBUG)
        logger.info("download %s", "ready")
        assert len(records) == 1
        assert records[0].getMessage() == "download ready"
        configured = next(item for item in logger.handlers if isinstance(item, logging.StreamHandler))
        assert configured.formatter.format(records[0]) == "[libephemeris] INFO: download ready"
        config.disable_logging()
        logger.critical("suppressed")
        assert len(records) == 1
    finally:
        logger.removeHandler(handler)


@pytest.mark.parametrize("size, expected", [
    (0, "0.0 B"), (512, "512.0 B"), (1023, "1023.0 B"),
    (1024, "1.0 KB"), (1536, "1.5 KB"),
    (119537664, "114.0 MB"), (1024**5, "1.0 PB"),
    (1024**6, "1024.0 PB"), (-5, "-5.0 B"), (1.25, "1.2 B"),
])
def test_file_size_binary_units(size, expected):
    assert config.format_file_size(size) == expected


def test_standard_scalar_errors():
    with pytest.raises(TypeError):
        config.set_log_level(object())
    with pytest.raises(TypeError):
        config.format_file_size(object())
