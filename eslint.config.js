// ESLint configuration (development only: skill users never need this).
// Run: npm run lint   ·   npm run lint:fix
import js from '@eslint/js';
import globals from 'globals';
import html from 'eslint-plugin-html';

const quality = {
  // correctness
  'no-unused-vars': ['error', { args: 'after-used', argsIgnorePattern: '^_', caughtErrors: 'none' }],
  'no-shadow': ['error', { builtinGlobals: false, hoist: 'functions' }],
  'no-use-before-define': ['error', { functions: false, classes: true, variables: false }], // closures may reference later consts
  'no-implicit-coercion': ['error', { allow: ['!!', '+'] }],
  'no-param-reassign': ['error', { props: false }],
  'no-return-assign': ['error', 'except-parens'],
  'no-self-compare': 'error',
  'no-template-curly-in-string': 'error',
  'no-unmodified-loop-condition': 'error',
  'no-unreachable-loop': 'error',
  'no-constant-binary-expression': 'error',
  'array-callback-return': 'error',
  eqeqeq: ['error', 'always', { null: 'ignore' }],
  // style that prevents bugs
  'no-var': 'error',
  'prefer-const': ['error', { destructuring: 'all' }],
  'object-shorthand': ['error', 'properties'],
  'prefer-template': 'error',
  'no-useless-concat': 'error',
  'no-else-return': ['error', { allowElseIf: false }],
  'no-lonely-if': 'error',
  curly: ['error', 'multi-line'],
  'dot-notation': 'error',
  'no-nested-ternary': 'off',          // dense numeric code reads fine with simple nested ternaries
};

export default [
  { ignores: ['node_modules/**', 'examples/**/*.html', '!examples/**/page.html', '**/*-shots/**'] },
  js.configs.recommended,
  {
    files: ['**/*.mjs', '**/*.js'],
    languageOptions: { ecmaVersion: 2024, sourceType: 'module', globals: { ...globals.browser, ...globals.node } },
    rules: quality,
  },
  {
    // Page sources (templates and example pages): lint their inline module scripts.
    files: ['templates/**/*.html', 'examples/**/page.html'],
    plugins: { html },
    languageOptions: { ecmaVersion: 2024, sourceType: 'module', globals: { ...globals.browser } },
    rules: quality,
  },
];
