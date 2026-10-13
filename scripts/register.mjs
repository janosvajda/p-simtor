// Usage: node --import <skill>/scripts/register.mjs --test physics.test.mjs
// Lets project files import "psim:core/..." and "psim:ui/..." modules from the skill.
import { register } from 'node:module';

register('./lib/psim-hooks.mjs', import.meta.url);
