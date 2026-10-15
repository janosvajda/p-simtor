import { test } from 'node:test';
import assert from 'node:assert/strict';
import { julianDay, sunEquatorial, sunPosition, sunDirection, sunriseSunset } from '../core/sun.mjs';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} ${a} vs ${b} (tol ${tol})`);
const hm = d => d.getUTCHours() * 60 + d.getUTCMinutes() + d.getUTCSeconds() / 60;

test('Julian day and Meeus example 25.a (1992 Oct 13, 0h): δ = −7.785°, α = 198.381°', () => {
  close(julianDay(new Date(Date.UTC(2000, 0, 1, 12))), 2451545, 1e-9, 'J2000');
  const s = sunEquatorial(new Date(Date.UTC(1992, 9, 13)));
  close(s.declination, -7.78507, 0.002, 'declination (Meeus low-precision result)');
  close(s.rightAscension, 198.38083, 0.002, 'right ascension');
});

test('equinox and solstice declinations; equation-of-time extremes', () => {
  close(sunEquatorial(new Date(Date.UTC(2024, 2, 20, 3, 6))).declination, 0, 0.01, 'March equinox 2024 (USNO 03:06 UTC)');
  close(sunEquatorial(new Date(Date.UTC(2024, 5, 20, 20, 51))).declination, 23.44, 0.01, 'June solstice 2024');
  close(sunEquatorial(new Date(Date.UTC(2024, 10, 3, 12))).equationOfTime, 16.4, 0.2, 'early November maximum');
  close(sunEquatorial(new Date(Date.UTC(2024, 1, 11, 12))).equationOfTime, -14.2, 0.2, 'mid-February minimum');
});

test('sky position: noon altitude 90° − φ + δ, the Sun due south at local solar noon, unit direction', () => {
  const date = new Date(Date.UTC(2024, 5, 21)), lat = 51.5, lon = -0.13, { solarNoon } = sunriseSunset(date, lat, lon);
  const p = sunPosition(solarNoon, lat, lon), d = sunEquatorial(solarNoon).declination;
  close(p.altitude, 90 - lat + d, 0.01, 'noon altitude'); close(p.azimuth, 180, 0.05, 'due south'); close(p.hourAngle, 0, 0.05);
  close(Math.hypot(...sunDirection(solarNoon, lat, lon)), 1, 1e-12);
});

test('sunrise and sunset match published times; polar day and night', () => {
  // London, 21 June 2024 (timeanddate.com): sunrise 04:43 BST, sunset 21:21 BST
  const ldn = sunriseSunset(new Date(Date.UTC(2024, 5, 21)), 51.5074, -0.1278);
  close(hm(ldn.sunrise), 3 * 60 + 43, 2, 'sunrise UTC'); close(hm(ldn.sunset), 20 * 60 + 21, 2, 'sunset UTC');
  close(ldn.dayLength, 16 + 38 / 60, 0.05, 'day length 16 h 38 min');
  assert.equal(sunriseSunset(new Date(Date.UTC(2024, 5, 21)), 78.2, 15.6).polar, 'day', 'Svalbard midsummer');
  assert.equal(sunriseSunset(new Date(Date.UTC(2024, 11, 21)), 78.2, 15.6).polar, 'night', 'Svalbard midwinter');
  close(sunriseSunset(new Date(Date.UTC(2024, 2, 20)), 0, 0).dayLength, 12.12, 0.03, 'equator: slightly over 12 h (refraction + disc)');
});
