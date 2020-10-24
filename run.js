
const PHRASE = 'beef';
const CMD = `./profanity.x64 --contract --skip 1 --matching ${PHRASE} -I 1000 -i 255`;

const exec = require('child_process').exec;
const fs = require('fs');

const loop = () => {
  const proc = exec(CMD);

  const { toChecksumAddress } = require('ethereum-checksum-address');
  const beep = require('beepbeep');

  const check = (priv, addr) => {
    console.log(`check: ${addr}`);
    if (addr.indexOf(PHRASE) !== -1) {
      fs.appendFileSync('output.txt', [toChecksumAddress(addr), priv].join('\t') + '\n');
      beep();
      proc.kill();
      loop();
    }
  };

  const match = str => {
    console.log(str);
    const m = str.match(/Private: (0x[a-f0-9]+) Contract: (0x[a-f0-9]+)/);
    if (m) check(m[1], m[2]);
  };

  proc.stdout.on('data', data => match(data + ''));
}

loop();
