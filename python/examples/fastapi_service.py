"""FastAPI example for Kerykeion Python 6 and the installed native v4 wheel.

Set KERYKEION_BACKEND=libephemeris, KERYKEION_LEB_MODE=leb and the native
data directory before launch. Run one server process, for example:
    uvicorn fastapi_service:app --app-dir python/examples --workers 1

FastAPI/Uvicorn/Kerykeion are application dependencies, not engine dependencies.
The existing Kerykeion session lock remains responsible for its global setters.
"""
from __future__ import annotations

from contextlib import asynccontextmanager
import os
import time
from typing import Literal

from fastapi import FastAPI, HTTPException
from fastapi.responses import JSONResponse, Response
from pydantic import BaseModel

from bounded_workers import BoundedWorkers, ServerBusyError


class SubjectRequest(BaseModel):
    """Small offline request used by the example and its load diagnostic."""

    year: int = 1990
    month: int = 6
    day: int = 15
    hour: int = 12
    minute: int = 0
    lat: float = 41.9028
    lng: float = 12.4964
    tz_str: str = "Europe/Rome"
    zodiac_type: Literal["Tropical", "Sidereal"] = "Tropical"
    perspective_type: Literal["Apparent Geocentric", "Topocentric"] = "Apparent Geocentric"
    default_points: bool = False


def initialize():
    """Initialize local native data before readiness; refuse the wrong backend."""
    from importlib.metadata import version
    from kerykeion.ephemeris_backend import BACKEND_NAME, ephe

    if BACKEND_NAME != "libephemeris" or not version("kerykeion").startswith("6."):
        raise RuntimeError("this example requires Kerykeion Python 6 with LibEphemeris")
    ephe.set_network_policy("sealed")
    ephe.set_calc_mode("leb")
    ephe.set_iers_delta_t_enabled(False)
    ephe.set_auto_spk_download(False)
    ephe.set_iers_auto_download(False)
    build_response(SubjectRequest(), False)


def build_response(request: SubjectRequest, svg: bool) -> Response:
    """Build the whole subject and encoded response on one calculating thread."""
    from kerykeion import AstrologicalSubjectFactory, ChartDataFactory, ChartDrawer

    keywords = request.model_dump(exclude={"default_points"})
    keywords.update(online=False, suppress_geonames_warning=True)
    if request.zodiac_type == "Sidereal":
        keywords["sidereal_mode"] = "LAHIRI"
    if not request.default_points:
        keywords["active_points"] = [
            "Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn",
            "Uranus", "Neptune", "Pluto", "Ascendant", "Medium_Coeli",
        ]
    subject = AstrologicalSubjectFactory.from_birth_data("FastAPI example", **keywords)
    if svg:
        chart = ChartDataFactory.create_chart_data("Natal", subject)
        return Response(ChartDrawer(chart, theme="classic").generate_svg_string(),
                        media_type="image/svg+xml")
    return JSONResponse(subject.model_dump(mode="json"))


@asynccontextmanager
async def lifespan(app: FastAPI):
    workers = int(os.environ.get("EPHEMERIS_WORKERS", "1"))
    queued = int(os.environ.get("EPHEMERIS_MAX_QUEUED", str(2 * workers)))
    async with BoundedWorkers(workers, queued) as pool:
        await pool.run(initialize)
        app.state.calculations = pool
        yield


app = FastAPI(lifespan=lifespan)


@app.get("/health")
async def health():
    return {"ready": True}


async def calculate(request: SubjectRequest, svg: bool):
    admitted = time.perf_counter_ns()

    def job():
        started = time.perf_counter_ns()
        response = build_response(request, svg)
        response.headers["X-Queue-Ns"] = str(started - admitted)
        response.headers["X-Work-Ns"] = str(time.perf_counter_ns() - started)
        return response

    try:
        return await app.state.calculations.run(job)
    except ServerBusyError as error:
        raise HTTPException(503, "calculation capacity is full", headers={"Retry-After": "1"}) from error


@app.post("/subject")
async def subject(request: SubjectRequest):
    return await calculate(request, False)


@app.post("/pipeline")
async def pipeline(request: SubjectRequest):
    return await calculate(request, True)
