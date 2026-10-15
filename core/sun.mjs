// Where the Sun is in the sky, for any date and place: the low-precision solar theory of Meeus,
// Astronomical Algorithms (2nd ed.), ch. 25, as used by NOAA's solar calculator. Accuracy about
// 0.01° in position and a few seconds in the equation of time for 1950–2050. Angles in degrees;
// dates are JavaScript Date objects (UTC). Atmospheric refraction is not applied except in the
// conventional sunrise/sunset altitude h₀ = −0.833° (refraction plus the Sun's radius).
//
//   julianDay(date)                      → JD (UT)
//   sunEquatorial(date)                  → { declination, rightAscension, eclipticLongitude, obliquity, equationOfTime (min) }
//   sunPosition(date, lat, lon)          → { altitude, azimuth (from north, clockwise), hourAngle }
//   sunDirection(date, lat, lon)         → unit vector [east, north, up]
//   sunriseSunset(date, lat, lon, h0)    → { sunrise, sunset, solarNoon (Date), dayLength (h), polar: 'day'|'night'|null }

const rad = Math.PI / 180, deg = 180 / Math.PI;
const norm360 = x => ((x % 360) + 360) % 360;

export function julianDay(date) {
  return date.getTime() / 86400000 + 2440587.5;
}

export function sunEquatorial(date) {
  const T = (julianDay(date) - 2451545) / 36525;
  const L0 = norm360(280.46646 + T * (36000.76983 + 0.0003032 * T));
  const M = 357.52911 + T * (35999.05029 - 0.0001537 * T);
  const e = 0.016708634 - T * (0.000042037 + 0.0000001267 * T);
  const C = (1.914602 - T * (0.004817 + 0.000014 * T)) * Math.sin(M * rad) + (0.019993 - 0.000101 * T) * Math.sin(2 * M * rad) + 0.000289 * Math.sin(3 * M * rad);
  const omega = 125.04 - 1934.136 * T;
  const lambda = L0 + C - 0.00569 - 0.00478 * Math.sin(omega * rad);
  const eps0 = 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - T * 0.001813))) / 60) / 60;
  const eps = eps0 + 0.00256 * Math.cos(omega * rad);
  const declination = Math.asin(Math.sin(eps * rad) * Math.sin(lambda * rad)) * deg;
  const rightAscension = norm360(Math.atan2(Math.cos(eps * rad) * Math.sin(lambda * rad), Math.cos(lambda * rad)) * deg);
  const y = Math.tan(eps * rad / 2) ** 2;
  const E = y * Math.sin(2 * L0 * rad) - 2 * e * Math.sin(M * rad) + 4 * e * y * Math.sin(M * rad) * Math.cos(2 * L0 * rad)
    - 0.5 * y * y * Math.sin(4 * L0 * rad) - 1.25 * e * e * Math.sin(2 * M * rad);
  return { declination, rightAscension, eclipticLongitude: norm360(lambda), obliquity: eps, equationOfTime: 4 * E * deg };
}

/** Altitude and azimuth (from north through east) of the Sun's centre, without refraction. */
export function sunPosition(date, lat, lon) {
  const { declination: d, equationOfTime } = sunEquatorial(date);
  const minutes = date.getUTCHours() * 60 + date.getUTCMinutes() + date.getUTCSeconds() / 60 + date.getUTCMilliseconds() / 60000;
  const hourAngle = (minutes + equationOfTime + 4 * lon) / 4 - 180;          // true solar time → hour angle
  const H = hourAngle * rad, phi = lat * rad, dr = d * rad;
  const altitude = Math.asin(Math.sin(phi) * Math.sin(dr) + Math.cos(phi) * Math.cos(dr) * Math.cos(H)) * deg;
  const azimuth = norm360(Math.atan2(Math.sin(H), Math.cos(H) * Math.sin(phi) - Math.tan(dr) * Math.cos(phi)) * deg + 180);
  return { altitude, azimuth, hourAngle: ((hourAngle + 540) % 360) - 180 };
}

/** Unit vector towards the Sun in local [east, north, up] coordinates (for shadows and lighting). */
export function sunDirection(date, lat, lon) {
  const { altitude, azimuth } = sunPosition(date, lat, lon);
  const a = altitude * rad, z = azimuth * rad;
  return [Math.cos(a) * Math.sin(z), Math.cos(a) * Math.cos(z), Math.sin(a)];
}

/**
 * Sunrise, sunset and solar noon (UTC) on the calendar day of `date` (UTC) at latitude/longitude,
 * for the altitude h0 (default −0.833°: upper limb on a refracted horizon). Iterates once with the
 * Sun's position at each event. polar is 'day' or 'night' when the Sun never sets or never rises.
 */
export function sunriseSunset(date, lat, lon, h0 = -0.833) {
  const day = Date.UTC(date.getUTCFullYear(), date.getUTCMonth(), date.getUTCDate());
  const at = minutes => new Date(day + minutes * 60000);
  const noonMin = m => 720 - 4 * lon - sunEquatorial(at(m)).equationOfTime;
  const solarNoon = noonMin(noonMin(720));
  const event = (sign, guess) => {
    let m = guess;
    for (let k = 0; k < 3; k++) {
      const d = sunEquatorial(at(m)).declination * rad, phi = lat * rad;
      const cosH = (Math.sin(h0 * rad) - Math.sin(phi) * Math.sin(d)) / (Math.cos(phi) * Math.cos(d));
      if (cosH < -1) return 'day';
      if (cosH > 1) return 'night';
      m = noonMin(m) + sign * 4 * Math.acos(cosH) * deg;
    }
    return m;
  };
  const rise = event(-1, solarNoon - 360), set = event(+1, solarNoon + 360);
  if (typeof rise === 'string' || typeof set === 'string') {
    const polar = rise === 'day' || set === 'day' ? 'day' : 'night';
    return { sunrise: null, sunset: null, solarNoon: at(solarNoon), dayLength: polar === 'day' ? 24 : 0, polar };
  }
  return { sunrise: at(rise), sunset: at(set), solarNoon: at(solarNoon), dayLength: (set - rise) / 60, polar: null };
}
