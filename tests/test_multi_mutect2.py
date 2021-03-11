#!/usr/bin/env python3

import pathlib
import subprocess
import unittest
from types import SimpleNamespace
from unittest import mock

from mutect2_tool import multi_mutect2 as MOD


class ThisTestCase(unittest.TestCase):
    def setUp(self):
        super().setUp()
        self.mocks = SimpleNamespace(
            subprocess=mock.MagicMock(spec_set=subprocess),
            futures=mock.MagicMock(spec_set=MOD.concurrent.futures),
        )

    def setup_popen(self, stdout=None, stderr=None, returncode=0):
        stdout = stdout or b''
        stderr = stderr or b''
        popen_instance = mock.create_autospec(subprocess.Popen, spec_set=True)
        self.mocks.subprocess.Popen.return_value = popen_instance
        popen_instance.communicate.return_value = (stdout, stderr)
        type(popen_instance).returncode = mock.PropertyMock(return_value=returncode)
        return popen_instance


class Test_subprocess_commands_pipe(ThisTestCase):
    def setUp(self):
        super().setUp()

    def tearDown(self):
        super().tearDown()

    def test_Popen_called_as_expected(self):
        mock_popen = self.setup_popen()
        cmd = MOD.CMD_STR
        MOD.subprocess_commands_pipe(cmd, di=self.mocks)
        self.mocks.subprocess.Popen.assert_called_once_with(
            MOD.shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        mock_popen.communicate.assert_called_once_with(timeout=3600)

    def test_ValueError_raised_with_TimeoutException(self):
        mock_popen = self.setup_popen()
        cmd = MOD.CMD_STR
        mock_popen.communicate.side_effect = [
            subprocess.TimeoutExpired(cmd, 3600),
            (b'', b'Timeout expired'),
        ]
        with self.assertRaises(ValueError):
            MOD.subprocess_commands_pipe(cmd, di=self.mocks)
        self.mocks.subprocess.Popen.assert_called_once_with(
            MOD.shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        expected_calls = [
            mock.call(timeout=3600),
            mock.call(),
        ]
        mock_popen.communicate.assert_has_calls(expected_calls)


class Test_ThreadPoolExecutor(ThisTestCase):
    def setUp(self):
        return super().setUp()

    def tearDown(self):
        return super().tearDown()

    def yield_futures(self, cmds):
        for cmd in cmds:
            yield mock.MagicMock(spec_set=MOD.concurrent.futures.Future)

    def test_commands_submitted(self):
        tpe_mock = mock.MagicMock(spec_set=MOD.concurrent.futures.ThreadPoolExecutor)
        self.mocks.futures.ThreadPoolExecutor.return_value.__enter__.return_value = (
            tpe_mock
        )
        commands = ['foo bar']
        mock_fn = mock.Mock()
        max_workers = 2
        MOD.tpe_submit_commands(
            commands, thread_count=max_workers, fn=mock_fn, di=self.mocks
        )
        expected_calls = [mock.call(mock_fn, cmd) for cmd in commands]
        tpe_mock.submit.assert_has_calls(expected_calls)
        self.mocks.futures.ThreadPoolExecutor.assert_called_once_with(
            max_workers=max_workers
        )

    @unittest.skip("Skipping for testing purposes")
    def test_submit_returns_futures(self):
        commands = list('abcde')
        tpe_mock = mock.MagicMock(spec_set=MOD.concurrent.futures.ThreadPoolExecutor)
        cmd_futures = self.yield_futures(commands)
        self.mocks.futures.as_completed.return_value = cmd_futures
        self.mocks.futures.ThreadPoolExecutor.return_value.__enter__.return_value = (
            tpe_mock
        )
        mock_fn = mock.Mock()
        max_workers = 2
        MOD.tpe_submit_commands(
            commands, thread_count=max_workers, fn=mock_fn, di=self.mocks
        )
        expected_calls = [mock.call(mock_fn, cmd) for cmd in commands]
        tpe_mock.submit.assert_has_calls(expected_calls)
        for future in cmd_futures:
            with self.subTest(future=future):
                future.result.assert_called_once_with()


# __END__
