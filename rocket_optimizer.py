#!/usr/bin/env python3

# TODO: Maybe have a Propulsion class that has a list of parts and is
# separate from SRB/LiquidEngine?  Things are getting hairy with
# decouplers adding decouple_direction field to SRBs/LiquidEngines.
# Methods of SRB/LiquidEngine need to preserve it, which isn't obvious
# from the definition of the SRB/LiquidEngine classes.

import argparse
import math
from enum import Enum
from copy import copy

ACCELERATION_GRAVITY = 9.81

parser = argparse.ArgumentParser(description='Help for finding the most efficient way to get delta-V!')
parser.add_argument('payload_mass', metavar='M', type=float,
                   help='Mass in tonnes of payload that these stages will deliver.')

parser.add_argument('--steerable1st', dest='steerable1st', action='store_true')
parser.add_argument('--no-steerable1st', dest='steerable1st', action='store_false')
parser.set_defaults(steerable1st=False)

parser.add_argument('--delta-v-of-atmosphere', dest='delta_v_of_atmosphere', type=float,
                    help='All stages that start below this delta_v will use atmosphere Isp.  Other will use vaccum.', default=1000)

parser.add_argument('--filter', dest='filter', action='store_true')
parser.add_argument('--no-filter', dest='filter', action='store_false')
parser.set_defaults(filter=True)

parser.add_argument('--min-twr', dest='min_twr_at_launch', type=float,
                    help='Minimum Thrust To Weight Ratio, i.e. acceleration as a multiple of g, at launch.',
                    default=1.2)

args = parser.parse_args()

class Radius(Enum):
    tiny = 1
    small = 2
    large = 3
    extra_large = 4
    radial = 5

class Direction(Enum):
    stack = 1
    radial = 2

class Part:
    def __init__(self, name, radius, cost, mass):
        self.name = name
        self.radius = radius
        self.cost = cost
        self.mass = mass
    def __radd__(self, other):
        c = copy(other)
        iadd(c, self)
        return c

# Also for separators
class Decoupler(Part):
    pass

def iadd(self, other):
    self.cost += other.cost
    self.full_name += " w/ " + other.name

    if hasattr(self, 'mass'):
        self.mass += other.mass
    else:
        self.full_mass += other.mass
        self.empty_mass += other.mass

    if isinstance(other, Decoupler):
        dir = Direction.radial if other.radius == Radius.radial else Direction.stack
        assert not hasattr(self, 'decouple_direction') or self.decouple_direction == dir
        self.decouple_direction = dir

SmallStackDecoupler = Decoupler("Small Stack Decoupler", Radius.small, 400, 0.05)
RadialDecoupler = Decoupler("TT-38K Radial Decoupler", Radius.radial, 600, 0.025)
# Unfortunatley, the game has a part called "Small Nose Cone," whose
# radius is tiny, while the "Aerodynamic Nose Cone"'s radius is small
AerodynamicNoseCone = Part("Aerodynamic Nose Cone", Radius.small, 240, 0.03)

class SRB:
    def __init__(self, name, cost, full_mass, empty_mass, thrust_atm, isp_atm, isp_vac):
        self.name = name
        self.full_name = name
        self.cost = cost
        self.full_mass = full_mass
        self.empty_mass = empty_mass
        self.thrust_atm = thrust_atm
        self.ve_atm = isp_atm * ACCELERATION_GRAVITY
        self.ve_vac = isp_vac * ACCELERATION_GRAVITY

    def __iadd__(self, part):
        assert isinstance(part, Part)
        assert part.radius == Raidus.radial or part.radius == Radius.small
        iadd(self, part)

    def __mul__(self, multiplier):
        c = copy(self)
        c.name += 'x' + str(multiplier)
        c.cost *= multiplier
        c.full_mass *= multiplier
        c.empty_mass *= multiplier
        c.thrust_atm *= multiplier
        # Don't need to modify ve_atm, ve_vac or decouple_direction
        return c

Flea = SRB('Flea', 200, 1.5, 0.45, 162.91, 140, 165)
Hammer = SRB('Hammer', 400, 3.56, 0.75, 197.90, 170, 195)
Thumper = SRB('Thumper', 850, 7.65, 1.5, 250.0, 175, 210)
Kickback = SRB('Kickback', 2700, 24, 4.5, 593.86, 195, 220)

FleaS = Flea + SmallStackDecoupler
HammerS = Hammer + SmallStackDecoupler
ThumperS = Thumper + SmallStackDecoupler
KickbackS = Kickback + SmallStackDecoupler
Thumperx3S = Thumper * 3 + AerodynamicNoseCone + AerodynamicNoseCone + SmallStackDecoupler

srbs = [FleaS, HammerS, ThumperS, KickbackS, Thumperx3S]

def fuel_cost(fuel, radius):
    assert radius is Radius.small or radius is Radius.large

    if radius is Radius.small:
        assert fuel % 100 == 0

        num800 = fuel // 800
        cost = 800 * num800
        fuel -= 800 * num800

        if fuel >= 400:
            cost += 500
            fuel -= 400

        if fuel >= 200:
            cost += 275
            fuel -= 200

        if fuel >= 100:
            cost += 150
            fuel -= 100

        assert fuel == 0

    elif radius is Radius.large:
        assert fuel % 800 == 0

        num6400 = fuel // 6400
        cost = 5750 * num6400
        fuel -= 6400 * num6400

        if fuel >= 3200:
            cost += 3000
            fuel -= 3200

        if fuel >= 1600:
            cost += 1550
            fuel -= 1600

        if fuel >= 800:
            cost += 800
            fuel -= 800

        assert fuel == 0

    return cost

class LiquidEngine:
    def __init__(self, name, radius, cost, mass, thrust_atm, isp_atm, isp_vac):
        self.name = name
        self.radius = radius
        self.cost = cost
        self.mass = mass
        self.thrust_atm = thrust_atm
        self.ve_atm = isp_atm * ACCELERATION_GRAVITY
        self.ve_vac = isp_vac * ACCELERATION_GRAVITY

Terrier = LiquidEngine('Terrier', Radius.small, 390, 0.5, 14.73, 85, 345)
Swivel = LiquidEngine('Swivel', Radius.small, 1200, 1.5, 168.75, 270, 320)
Vector = LiquidEngine('Vector', Radius.small, 18000, 4.0, 936.5, 295, 315)

Poodle = LiquidEngine('Poodle', Radius.large, 1300, 1.75, 64.29, 90, 350)
Skipper = LiquidEngine('Skipper', Radius.large, 5300, 3.0, 568.75, 280, 320)
Mainsail = LiquidEngine('Mainsail', Radius.large, 13000, 6.0, 1379.0, 285, 310)

# Liquid engine plus fuel tanks, plus optional SRBs.
class Liquid:
    def __init__(self, engine, fuel, srbs = None):
        if srbs is None:
            srbs = []

        self.engine = engine
        self.fuel = fuel
        # For now, all liquid stages use stack decouplers.
        self.decouple_direction = Direction.stack

        # Compute weighted average of exhaust velocities.  I think
        # this is accurate if we throttle everything such that they
        # all burn out at the same time.

        # For liquid fuel tanks, the weight (both full and empty) is
        # proportional to capacity.
        liquid_fuel_mass = fuel / 200.0
        ve_atm = engine.ve_atm * liquid_fuel_mass
        ve_vac = engine.ve_vac * liquid_fuel_mass
        total_fuel_mass = liquid_fuel_mass
        for booster in srbs:
            booster_fuel_mass = booster.full_mass - booster.empty_mass
            ve_atm += booster.ve_atm * booster_fuel_mass
            ve_vac += booster.ve_vac * booster_fuel_mass
            total_fuel_mass += booster_fuel_mass

        self.ve_atm = ve_atm / total_fuel_mass
        self.ve_vac = ve_vac / total_fuel_mass

        self.name = engine.name + " " + str(fuel)
        for booster in srbs:
            self.name += " " + booster.name

        liquid_empty_mass = engine.mass + fuel / 1600.0
        self.empty_mass = liquid_empty_mass + sum([s.empty_mass for s in srbs])
        self.full_mass = liquid_empty_mass + liquid_fuel_mass + \
            sum([s.full_mass for s in srbs])
        self.thrust_atm = engine.thrust_atm + sum([s.thrust_atm for s in srbs])
        self.cost = engine.cost + fuel_cost(fuel, engine.radius) + \
            sum([s.cost for s in srbs])


class Stage:
    def __init__(self, rocket, propulsion, payload_mass):
        assert hasattr(propulsion, 'decouple_direction')
        self.rocket = rocket
        self.propulsion = propulsion
        self.cost = propulsion.cost
        self.payload_mass = payload_mass
        self.full_mass = payload_mass + propulsion.full_mass
        self.empty_mass = payload_mass + propulsion.empty_mass

    def compute_delta_v(self, atm):
        self.ve = self.propulsion.ve_atm if atm else self.propulsion.ve_vac
        self.delta_v = self.ve * math.log(self.full_mass / self.empty_mass)

class Rocket:
    # Stages ordered from highest to lowest, i.e. from closest to
    # payload to furthest from payload.
    def __init__(self, propulsions):
        propulsions = [p for p in propulsions if p is not None]

        # Loop through stages from payload to ground, computing mass of each stage as we go.
        self.stages = []
        prev_mass = args.payload_mass
        for i, propulsion in enumerate(propulsions):
            stage = Stage(self, propulsion, prev_mass)
            self.stages.append(stage)
            prev_mass = stage.full_mass

        # Can compute TWR at launch
        thrust = 0
        for stage in reversed(self.stages):
            thrust += stage.propulsion.thrust_atm
            if stage.propulsion.decouple_direction == Direction.stack:
                break

        self.twr_launch = thrust / self.stages[-1].full_mass / ACCELERATION_GRAVITY

        # Now loop in the other direction, deciding when to use atm of
        # vac Ve, and computing delta_v
        prev_delta_v = 0
        for s in reversed(self.stages):
            s.compute_delta_v(prev_delta_v < args.delta_v_of_atmosphere)
            prev_delta_v += s.delta_v

        self.cost = sum([s.cost for s in self.stages])
        self.delta_v = sum([s.delta_v for s in self.stages])

# Up to two stages, all SRBs:
#rockets = \
#    [Rocket([first]) for first in srbs] + \
#    [Rocket([second, first]) for second in srbs for first in srbs] + \
#    [Rocket([third, second, first]) for third in srbs for second in srbs for first in srbs]

Thumper_radial = Thumper + AerodynamicNoseCone
ThumperR = Thumper_radial + RadialDecoupler

terriers = [Liquid(Terrier, fuel) for fuel in range(100, 3201, 100)]
swivels = [Liquid(Swivel, fuel) for fuel in range(100, 3201, 100)]
swivels_srbs = [Liquid(Swivel, fuel, [Thumper_radial, Thumper_radial]) for fuel in range(100, 3201, 100)]

radial_srbs = [ThumperR * 2]

#terriers = [Liquid(Terrier, 1600)]
#swivels = [Liquid(Swivel, 1600)]
#swivels_srbs = [Liquid(Swivel, 1600, [Thumper_radial, Thumper_radial])]

# Swivel is never chosen for the middle of three stages.
#
# For a 4.53t payload, it never makes sense to "split" the terrier
# stage into two stages, since the cost of the terrier + stack
# separator = 790, about the cost of an FL-T800 (800).  So losing 800
# fuel is just not worth it.
#
# For a 0.53t payload, it starts to make sense around 6700 m/s delta
# v (atmosphere at 500 m/s).  But who needs that much delta v??

rockets = \
    [Rocket([single]) for single in swivels] + \
    [Rocket([second, first])
         for second in terriers
         for first in swivels + swivels_srbs] + \
    [Rocket([third, second, first])
         for third in terriers
         for second in terriers + swivels
         for first in swivels + swivels_srbs] + \
    [Rocket([third, second, first])
         for third in terriers
         for second in swivels
         for first in radial_srbs] + \
    [Rocket([fourth, third, second, first])
         for fourth in terriers
         for third in terriers
         for second in swivels
         for first in radial_srbs]

# Get rid of ones that won't get off the launch pad fast enough.
rockets = [r for r in rockets if r.twr_launch >= args.min_twr_at_launch]

# Could add this just to eliminate silly ones.
# rockets = [r for r in rockets if r.stages[-1].delta_v >= args.min_delta_v_of_first_stage]

rockets.sort(key=lambda r : r.delta_v, reverse=True)
rockets.sort(key=lambda r : r.cost)

# Keep only the dominating ones, i.e. eliminate ones that cost the same or more, but have lower delta-v.

if args.filter:
    filtered_rockets = []
    best_delta_v = None
    for r in rockets:
        if best_delta_v is None or r.delta_v > best_delta_v:
            filtered_rockets.append(r)
            best_delta_v = r.delta_v
else:
    filtered_rockets = rockets

for r in filtered_rockets:
    print("%5d %6.1f %4.2f" % (r.cost, r.delta_v, r.twr_launch), 
          ["%10s, %7.2f, %7.2f" % (s.propulsion.name, s.delta_v, s.ve) for s in r.stages])
