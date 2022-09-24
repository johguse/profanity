
const PARSED = require('dotenv').config().parsed;
const {
  PHRASE
} = PARSED;
const CONTRACT = PARSED.CONTRACT === 'true';

const CMD = `./profanity.x64 ${CONTRACT ? '--contract' : ''} --skip 1 --matching ${PHRASE} -I 200 -i 255`;
console.log(CMD);

const exec = require('child_process').exec;
const fs = require('fs');

const loop = () => {
  exec('rm cache-opencl**');
  const proc = exec(CMD);

  const { toChecksumAddress } = require('ethereum-checksum-address');
  const beep = require('beepbeep');

  const check = (priv, addr) => {
    const checksum = toChecksumAddress(addr);
    console.log(`checksum: ${checksum}`);
    if (addr.indexOf(PHRASE.toLowerCase()) !== -1) {
      const some = [PHRASE, PHRASE.toLowerCase(), PHRASE.toUpperCase()].some(s => checksum.indexOf(s) !== -1);
      if (some) {
        const type = CONTRACT ? 'contract' : 'address';
        fs.appendFileSync('output.txt', [checksum, priv, type].join('\t') + '\n');
        beep();
      }
      beep();
      console.log('restart');
      proc.kill();
      setTimeout(loop, 1);
    }
  };

  const match = str => {
    console.log(str);
    const m = (CONTRACT ?
      str.match(/Private: (0x[a-f0-9]+) Contract: (0x[a-f0-9]+)/) :
      str.match(/Private: (0x[a-f0-9]+) Address: (0x[a-f0-9]+)/)
    );
    if (m) check(m[1], m[2]);
  };

  proc.stdout.on('data', data => match(data + ''));
}

loop();
